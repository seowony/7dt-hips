#!/usr/bin/env python3

import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.io import fits
from mocpy import MOC


INPUT_ROOT = Path("./hips")
OUTPUT_PATH = INPUT_ROOT / "color"
TILE_WIDTH = 512
OVERWRITE = True

# Leave END_ORDER to None to use the maximum order found in the input MOCs.
START_ORDER = 0
END_ORDER = None


# ----------------------------------------------------------------------
# Spectral weighting
# ----------------------------------------------------------------------

WEIGHTS = None

RED_WEIGHT_ROLLOFF_START_NM = 650
RED_WEIGHT_AT_MAX_WAVELENGTH = 1.00


# ----------------------------------------------------------------------
# Per-band pedestal
# ----------------------------------------------------------------------

PEDESTAL_FRACTION_AT_MIN_WAVELENGTH = 0.03
PEDESTAL_FRACTION_AT_MIN_WAVELENGTH = 0.04
PEDESTAL_RED_RAMP_START_NM = 625
PEDESTAL_FRACTION_AT_MAX_WAVELENGTH = 0.03
PEDESTAL_FRACTION_AT_MAX_WAVELENGTH = 0.04


# ----------------------------------------------------------------------
# RGB balance
# ----------------------------------------------------------------------

RGB_CHANNEL_GAINS = np.array(
#[1.02, 0.93, 0.97],
    [1.09, 0.90, 1.04],
    dtype=np.float32,
)
RGB_CHANNEL_GAINS /= 1.5

RGB_BLACKPOINT = np.array(
    [0.032, 0.010, 0.011],
    dtype=np.float32,
)

ORDER_BLACKPOINT_START = 8
ORDER_BLACKPOINT_FULL = 10

ORDER_BLACKPOINT_AT_FULL = np.array(
    [0.004, 0.004, 0.004],
    dtype=np.float32,
)


# ----------------------------------------------------------------------
# Stretch
# ----------------------------------------------------------------------

PER_BAND_NORM = "hips_pixel_cut_pedestal"
I_MODE = "luma"

Q = 7.0
Q = 4.0
Q = 15.0
I0_MODE = 4.5

GAMMA = 1.25
FINAL_DENOM = 8.5
FINAL_DENOM = 12


# ----------------------------------------------------------------------
# Saturation
# ----------------------------------------------------------------------
#
# Dark background remains relatively desaturated.
# Brighter sources get more colour.

SATURATION_SHADOW = 0.14
SATURATION_HIGHLIGHT = 0.80

SATURATION_LUMA_LOW = 0.08
SATURATION_LUMA_HIGH = 0.45


POOL_SIZE = 16


# ======================================================================
# Bands / properties
# ======================================================================

def discover_bands(root: Path) -> list[dict]:
    bands = []
    pattern = re.compile(r"^m(\d+)$")

    for path in root.iterdir():
        if not path.is_dir():
            continue

        match = pattern.match(path.name)
        if not match:
            continue

        moc_path = path / "Moc.fits"
        if not moc_path.exists():
            continue

        wavelength_nm = int(match.group(1))
        properties = read_properties(path / "properties")

        bands.append(
            {
                "id": path.name,
                "path": path,
                "moc_path": moc_path,
                "properties": properties,
                "wavelength_nm": wavelength_nm,
            }
        )

    bands.sort(key=lambda band: band["wavelength_nm"])

    if not bands:
        raise RuntimeError(f"No input HiPS found in {root}")

    return bands


def read_properties(path: Path) -> dict[str, str]:
    properties = {}

    if not path.exists():
        return properties

    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()

        if not line or line.startswith("#") or "=" not in line:
            continue

        key, value = line.split("=", 1)
        properties[key.strip()] = value.strip()

    return properties


def property_float_pair(
    properties: dict[str, str],
    key: str,
) -> tuple[float, float] | None:

    value = properties.get(key)

    if value is None:
        return None

    parts = value.split()

    if len(parts) < 2:
        return None

    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        return None


BANDS = discover_bands(INPUT_ROOT)

WAVELENGTHS_NM = np.array(
    [band["wavelength_nm"] for band in BANDS],
    dtype=np.float32,
)


# ======================================================================
# Weights / pedestal
# ======================================================================

def wavelength_ramp(
    values_nm: np.ndarray,
    start: float,
    end: float,
    start_value: float,
    end_value: float,
) -> np.ndarray:

    if end <= start:
        return np.full_like(
            values_nm,
            end_value,
            dtype=np.float32,
        )

    t = (values_nm - start) / (end - start)
    t = np.clip(t, 0.0, 1.0)

    return (
        start_value
        + t * (end_value - start_value)
    ).astype(np.float32)


def default_weights(
    wavelengths_nm: np.ndarray,
) -> np.ndarray:

    max_nm = float(np.max(wavelengths_nm))

    return wavelength_ramp(
        wavelengths_nm,
        RED_WEIGHT_ROLLOFF_START_NM,
        max_nm,
        1.0,
        RED_WEIGHT_AT_MAX_WAVELENGTH,
    )


def default_pedestal_fractions(
    wavelengths_nm: np.ndarray,
) -> np.ndarray:

    max_nm = float(np.max(wavelengths_nm))

    return wavelength_ramp(
        wavelengths_nm,
        PEDESTAL_RED_RAMP_START_NM,
        max_nm,
        PEDESTAL_FRACTION_AT_MIN_WAVELENGTH,
        PEDESTAL_FRACTION_AT_MAX_WAVELENGTH,
    )


WEIGHTS = (
    default_weights(WAVELENGTHS_NM)
    if WEIGHTS is None
    else np.asarray(WEIGHTS, dtype=np.float32)
)

if len(WEIGHTS) != len(BANDS):
    raise ValueError(
        "WEIGHTS must have one value per input HiPS band"
    )


BAND_PIXEL_CUTS = [
    property_float_pair(
        band["properties"],
        "hips_pixel_cut",
    )
    for band in BANDS
]

BAND_PEDESTAL_FRACTIONS = (
    default_pedestal_fractions(WAVELENGTHS_NM)
)


# ======================================================================
# RGB basis
# ======================================================================

def rgb_from_hue(h: float) -> np.ndarray:
    h = (h % 360.0) / 60.0
    c = 1.0
    x = 1.0 - abs(h % 2 - 1.0)

    if 0 <= h < 1:
        r, g, b = c, x, 0
    elif 1 <= h < 2:
        r, g, b = x, c, 0
    elif 2 <= h < 3:
        r, g, b = 0, c, x
    elif 3 <= h < 4:
        r, g, b = 0, x, c
    elif 4 <= h < 5:
        r, g, b = x, 0, c
    else:
        r, g, b = c, 0, x

    return np.array(
        [r, g, b],
        dtype=np.float32,
    )


def triangular(
    x: np.ndarray,
    center: float,
    half_width: float,
) -> np.ndarray:

    return np.clip(
        1.0 - np.abs(x - center) / half_width,
        0.0,
        1.0,
    )


def make_rgb_basis(
    wavelengths_nm: np.ndarray,
) -> np.ndarray:

    red = triangular(
        wavelengths_nm,
        760.0,
        190.0,
    )

    green = triangular(
        wavelengths_nm,
        575.0,
        155.0,
    )

    blue = triangular(
        wavelengths_nm,
        450.0,
        110.0,
    )

    basis = np.stack(
        [red, green, blue],
        axis=1,
    ).astype(np.float32)

    basis /= np.maximum(
        np.max(
            basis,
            axis=1,
            keepdims=True,
        ),
        1e-6,
    )

    # ONLY NEW CHANGE FOR THIS TEST:
    # equalize total contribution of R, G and B.
    channel_sum = np.sum(
        basis,
        axis=0,
        keepdims=True,
    )

    mean_sum = np.mean(channel_sum)

    basis *= mean_sum / np.maximum(
        channel_sum,
        1e-6,
    )

    return basis

RGB_BASIS = make_rgb_basis(WAVELENGTHS_NM)


# ======================================================================
# HiPS paths / FITS
# ======================================================================

def tile_suffix(
    order: int | None,
    ipix: int | None,
) -> str:

    if order is None and ipix is None:
        return "Norder3/Allsky.fits"

    return (
        f"Norder{order}/"
        f"Dir{int(10000 * (ipix // 10000))}/"
        f"Npix{ipix}.fits"
    )


def load_fits_any(
    path: Path,
) -> np.ndarray:

    """Return the first 2D FITS image found as float32."""

    with fits.open(
        path,
        memmap=True,
    ) as hdul:

        for key in (
            "SCI",
            "IMAGE",
            "PRIMARY",
        ):
            try:
                if (
                    key in hdul
                    and hdul[key].data is not None
                ):
                    return np.asarray(
                        hdul[key].data,
                        dtype=np.float32,
                    )

            except Exception:
                pass

        for hdu in hdul:
            data = getattr(
                hdu,
                "data",
                None,
            )

            if data is None:
                continue

            arr = np.asarray(data)

            if arr.ndim >= 2:
                return arr.astype(
                    np.float32
                )

    raise ValueError(
        f"No 2D image in {path}"
    )


def sanitize(
    x: np.ndarray,
) -> np.ndarray:

    x = np.asarray(
        x,
        dtype=np.float32,
    )

    if np.ma.isMaskedArray(x):
        x = x.filled(np.nan)

    return np.nan_to_num(
        x,
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )


# ======================================================================
# Per-band normalization
# ======================================================================

def normalize_per_band(
    bands: list[np.ndarray],
    method: str = "hips_pixel_cut_high",
    pixel_cuts: list[
        tuple[float, float] | None
    ] | None = None,
    pedestal_fractions: np.ndarray | None = None,
) -> list[np.ndarray]:

    out = []

    pixel_cuts = (
        pixel_cuts
        or [None] * len(bands)
    )

    if pedestal_fractions is None:
        pedestal_fractions = np.zeros(
            len(bands),
            dtype=np.float32,
        )

    for (
        band,
        pixel_cut,
        pedestal_fraction,
    ) in zip(
        bands,
        pixel_cuts,
        pedestal_fractions,
    ):

        if np.ma.isMaskedArray(band):
            valid = np.isfinite(
                band.filled(np.nan)
            )
        else:
            valid = np.isfinite(band)

        band = sanitize(band)

        if method == "none":

            denom = 1.0
            normalized = band

        elif method == "hips_pixel_cut_positive":

            if pixel_cut is None:
                raise ValueError(
                    "hips_pixel_cut_positive needs "
                    "hips_pixel_cut in every input properties file"
                )

            _, high = pixel_cut

            denom = high

            normalized = (
                np.clip(
                    band,
                    0.0,
                    high,
                )
                / denom
            )

        elif method == "hips_pixel_cut_pedestal":

            if pixel_cut is None:
                raise ValueError(
                    "hips_pixel_cut_pedestal needs "
                    "hips_pixel_cut in every input properties file"
                )

            _, high = pixel_cut

            pedestal = (
                pedestal_fraction
                * high
            )

            denom = (
                high
                - pedestal
            )

            normalized = (
                np.clip(
                    band - pedestal,
                    0.0,
                    denom,
                )
                / denom
            )

        elif method == "hips_pixel_cut_high":

            if pixel_cut is None:
                raise ValueError(
                    "hips_pixel_cut_high needs "
                    "hips_pixel_cut in every input properties file"
                )

            _, high = pixel_cut

            denom = high

            normalized = (
                band / denom
            )

        elif method == "hips_pixel_cut_range":

            if pixel_cut is None:
                raise ValueError(
                    "hips_pixel_cut_range needs "
                    "hips_pixel_cut in every input properties file"
                )

            low, high = pixel_cut

            denom = (
                high - low
            )

            normalized = (
                np.clip(
                    band,
                    low,
                    high,
                )
                - low
            ) / denom

        elif method == "meanpos":

            denom = (
                np.nanmean(
                    band[band > 0]
                )
                if np.any(band > 0)
                else 1.0
            )

            normalized = (
                band / denom
            )

        elif method == "max":

            denom = (
                np.nanmax(band)
                if np.size(band)
                else 1.0
            )

            normalized = (
                band / denom
            )

        elif method == "q999":

            denom = (
                np.nanquantile(
                    band[band > 0],
                    0.999,
                )
                if np.any(band > 0)
                else 1.0
            )

            normalized = (
                band / denom
            )

        else:

            raise ValueError(
                "Unknown per-band normalization "
                f"method: {method}"
            )

        if (
            not np.isfinite(denom)
            or denom <= 0
        ):
            denom = 1.0
            normalized = band

        normalized = np.where(
            valid,
            normalized,
            0.0,
        )

        out.append(
            normalized.astype(
                np.float32
            )
        )

    return out


# ======================================================================
# Composite
# ======================================================================

def composite_linear(
    bands: list[np.ndarray],
    weights: np.ndarray,
    rgb_basis: np.ndarray,
    per_band_norm: str = "meanpos",
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:

    normalized = normalize_per_band(
        bands,
        method=per_band_norm,
        pixel_cuts=BAND_PIXEL_CUTS,
        pedestal_fractions=BAND_PEDESTAL_FRACTIONS,
    )

    weighted = [
        weight * band
        for weight, band
        in zip(
            weights,
            normalized,
        )
    ]

    s_stack = np.stack(
        weighted,
        axis=-1,
    )

    rgb_lin = (
        s_stack
        @ rgb_basis
    )

    rgb_lin = np.clip(
        np.nan_to_num(
            rgb_lin,
            nan=0.0,
            posinf=0.0,
            neginf=0.0,
        ),
        0.0,
        None,
    )

    return (
        s_stack,
        rgb_lin,
        rgb_lin.copy(),
    )


def compute_intensity(
    s_stack: np.ndarray,
    rgb_lin: np.ndarray | None,
    mode: str = "sum",
    luma_weights=(
        0.2126,
        0.7152,
        0.0722,
    ),
) -> np.ndarray:

    if mode == "sum":

        intensity = np.clip(
            np.nansum(
                s_stack,
                axis=-1,
            ),
            0.0,
            None,
        )

    elif mode == "luma":

        if rgb_lin is None:
            raise ValueError(
                "rgb_lin is required for luma intensity"
            )

        rw, gw, bw = luma_weights

        intensity = np.clip(
            rw * rgb_lin[..., 0]
            + gw * rgb_lin[..., 1]
            + bw * rgb_lin[..., 2],
            0.0,
            None,
        )

    else:

        raise ValueError(
            f"Unknown intensity mode: {mode}"
        )

    return intensity


def lupton_scale(
    i: np.ndarray,
    q: float = 10.0,
    i0_mode: float = 1.0,
) -> np.ndarray:

    eps = 1e-6

    x = i / (
        i0_mode + eps
    )

    return (
        np.arcsinh(q * x)
        / np.maximum(
            x,
            eps,
        )
    ).astype(np.float32)


def composite_lupton(
    bands: list[np.ndarray],
    weights: np.ndarray,
    rgb_basis: np.ndarray,
    per_band_norm: str = "meanpos",
    I_mode: str = "sum",
    Q: float = 10.0,
    I0_mode: float = 1.0,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:

    (
        s_stack,
        rgb_lin,
        _,
    ) = composite_linear(
        bands,
        weights,
        rgb_basis,
        per_band_norm=per_band_norm,
    )

    intensity = compute_intensity(
        s_stack,
        rgb_lin,
        mode=I_mode,
    )

    scale = lupton_scale(
        intensity,
        q=Q,
        i0_mode=I0_mode,
    )

    rgb = (
        rgb_lin
        * scale[..., None]
    )

    rgb = np.clip(
        np.nan_to_num(
            rgb,
            nan=0.0,
            posinf=0.0,
            neginf=0.0,
        ),
        0.0,
        None,
    )

    return (
        s_stack,
        rgb_lin,
        rgb,
    )


# ======================================================================
# Display adjustments
# ======================================================================

def order_blackpoint(
    order: int | None,
) -> np.ndarray:

    if (
        order is None
        or ORDER_BLACKPOINT_FULL
        <= ORDER_BLACKPOINT_START
    ):
        return np.zeros(
            3,
            dtype=np.float32,
        )

    t = (
        order
        - ORDER_BLACKPOINT_START
    ) / (
        ORDER_BLACKPOINT_FULL
        - ORDER_BLACKPOINT_START
    )

    t = np.clip(
        t,
        0.0,
        1.0,
    )

    return (
        t
        * ORDER_BLACKPOINT_AT_FULL
    ).astype(np.float32)


def apply_rgb_channel_gains(
    rgb: np.ndarray,
    order: int | None,
) -> np.ndarray:

    rgb = (
        rgb
        * RGB_CHANNEL_GAINS
    )

    return np.clip(
        rgb
        - RGB_BLACKPOINT
        - order_blackpoint(order),
        0.0,
        None,
    )


def final_normalize_gamma_clip(
    rgb: np.ndarray,
    norm: str = "q999",
    q: float = 0.999,
    use_gamma: bool = True,
    gamma: float = 2.2,
    denom: float | None = None,
) -> np.ndarray:

    rgb = np.nan_to_num(
        rgb,
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )

    rgb = np.clip(
        rgb,
        0.0,
        None,
    )

    if denom is None:

        if norm == "q999":

            positive = rgb[
                rgb > 0
            ]

            denom = (
                np.nanquantile(
                    positive,
                    q,
                )
                if positive.size
                else 1.0
            )

        elif norm == "max":

            denom = np.nanmax(
                rgb
            )

        else:

            raise ValueError(
                "Unknown final normalization "
                f"method: {norm}"
            )

    if (
        not np.isfinite(denom)
        or denom <= 0
    ):
        denom = 1.0

    rgb = (
        rgb
        / (denom + 1e-6)
    )

    if use_gamma:

        rgb = np.power(
            rgb,
            1.0 / gamma,
            dtype=np.float32,
        )

    return np.clip(
        rgb,
        0.0,
        1.0,
    ).astype(np.float32)


def apply_saturation(
    rgb: np.ndarray,
) -> np.ndarray:

    """
    Luminance-dependent saturation.

    The sky background stays relatively neutral,
    while brighter sources retain more colour.
    """

    luma = (
        0.2126 * rgb[..., 0]
        + 0.7152 * rgb[..., 1]
        + 0.0722 * rgb[..., 2]
    )

    t = (
        luma
        - SATURATION_LUMA_LOW
    ) / (
        SATURATION_LUMA_HIGH
        - SATURATION_LUMA_LOW
    )

    t = np.clip(
        t,
        0.0,
        1.0,
    )

    # Smoothstep
    t = (
        t
        * t
        * (3.0 - 2.0 * t)
    )

    saturation = (
        SATURATION_SHADOW
        + t
        * (
            SATURATION_HIGHLIGHT
            - SATURATION_SHADOW
        )
    )

    out = (
        luma[..., None]
        + saturation[..., None]
        * (
            rgb
            - luma[..., None]
        )
    )

    return np.clip(
        out,
        0.0,
        1.0,
    ).astype(np.float32)


# ======================================================================
# PNG
# ======================================================================

def save_png(
    path: Path,
    rgb: np.ndarray,
) -> None:

    rgb_uint8 = (
        np.clip(
            rgb,
            0.0,
            1.0,
        )
        * 255.0
        + 0.5
    ).astype(np.uint8)

    try:
        import imageio.v2 as imageio

        rgb_uint8 = np.flipud(
            rgb_uint8
        )

        imageio.imwrite(
            path,
            rgb_uint8,
        )

    except Exception:
        from PIL import Image

        Image.fromarray(
            rgb_uint8,
            mode="RGB",
        ).save(path)


# ======================================================================
# Tile processing
# ======================================================================

def process(param):

    order, ipix = param

    suffix_path = tile_suffix(
        order,
        ipix,
    )

    out_tile_path = (
        OUTPUT_PATH
        / suffix_path.replace(
            ".fits",
            ".png",
        )
    )

    if (
        out_tile_path.exists()
        and not OVERWRITE
    ):
        return

    bands = []

    for band in BANDS:

        tile_path = (
            band["path"]
            / suffix_path
        )

        if tile_path.exists():

            bands.append(
                load_fits_any(
                    tile_path
                )
            )

        else:

            bands.append(
                np.full(
                    (
                        TILE_WIDTH,
                        TILE_WIDTH,
                    ),
                    np.nan,
                    dtype=np.float32,
                )
            )

    _, _, rgb = composite_lupton(
        bands,
        WEIGHTS,
        RGB_BASIS,
        per_band_norm=PER_BAND_NORM,
        I_mode=I_MODE,
        Q=Q,
        I0_mode=I0_MODE,
    )

    rgb = apply_rgb_channel_gains(
        rgb,
        order,
    )

    rgb = final_normalize_gamma_clip(
        rgb,
        use_gamma=True,
        gamma=GAMMA,
        denom=FINAL_DENOM,
    )

    # Saturation only affects RGB values.
    # It does not resample, blur or interpolate the image.
    rgb = apply_saturation(
        rgb
    )

    out_tile_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    save_png(
        out_tile_path,
        rgb,
    )

    print(
        f"Written: {out_tile_path}"
    )


# ======================================================================
# MOC / metadata
# ======================================================================

def build_union_moc() -> MOC:

    moc = None

    for band in BANDS:

        band_moc = MOC.from_fits(
            band["moc_path"]
        )

        moc = (
            band_moc
            if moc is None
            else moc.union(
                band_moc
            )
        )

    return moc


def has_input_allsky() -> bool:

    return any(
        (
            band["path"]
            / "Norder3/Allsky.fits"
        ).exists()
        for band in BANDS
    )


def write_output_metadata(
    moc: MOC,
    end_order: int,
) -> None:

    OUTPUT_PATH.mkdir(
        parents=True,
        exist_ok=True,
    )

    try:

        moc.save(
            str(
                OUTPUT_PATH
                / "Moc.fits"
            ),
            "fits",
            overwrite=True,
        )

    except TypeError:

        moc.save(
            str(
                OUTPUT_PATH
                / "Moc.fits"
            ),
            overwrite=True,
        )

    wavelengths = ", ".join(
        f'{band["id"]}='
        f'{band["wavelength_nm"]}nm'
        for band in BANDS
    )

    properties = [
        "creator_did          = ivo://UNK.AUTH/P/color",
        "obs_title            = 7DT color composite",

        (
            "obs_description      = "
            f"RGB composite from {len(BANDS)} "
            f"7DT HiPS bands: {wavelengths}"
        ),

        (
            "moc_sky_fraction     = "
            f"{moc.sky_fraction:.6E}"
        ),

        "hips_builder         = make_color_hips.py",
        "hips_version         = 1.4",
        "hips_frame           = equatorial",

        f"hips_order           =10",
        f"hips_order_min       = {START_ORDER}",
        f"hips_tile_width      = {TILE_WIDTH}",

        "hips_status          = public master clonableOnce",
        "hips_tile_format     = png",
        "dataproduct_type     = image",

        f"color_per_band_norm  = {PER_BAND_NORM}",
        f"color_i_mode         = {I_MODE}",
        f"color_q              = {Q}",
        f"color_i0             = {I0_MODE}",
        f"color_gamma          = {GAMMA}",
        f"color_final_denom    = {FINAL_DENOM}",

        (
            "color_saturation_shadow = "
            f"{SATURATION_SHADOW}"
        ),

        (
            "color_saturation_highlight = "
            f"{SATURATION_HIGHLIGHT}"
        ),

        (
            "color_saturation_luma_low = "
            f"{SATURATION_LUMA_LOW}"
        ),

        (
            "color_saturation_luma_high = "
            f"{SATURATION_LUMA_HIGH}"
        ),

        (
            "color_rgb_gains       = "
            + " ".join(
                f"{gain:g}"
                for gain
                in RGB_CHANNEL_GAINS
            )
        ),

        (
            "color_rgb_blackpoint  = "
            + " ".join(
                f"{black:g}"
                for black
                in RGB_BLACKPOINT
            )
        ),

        (
            "color_order_blackpoint= "
            f"start:{ORDER_BLACKPOINT_START} "
            f"full:{ORDER_BLACKPOINT_FULL} "
            + " ".join(
                f"{black:g}"
                for black
                in ORDER_BLACKPOINT_AT_FULL
            )
        ),

        (
            "color_red_weight      = "
            f"start_nm:{RED_WEIGHT_ROLLOFF_START_NM} "
            f"max_weight:{RED_WEIGHT_AT_MAX_WAVELENGTH}"
        ),

        (
            "color_pedestal_frac   = "
            f"start_nm:{PEDESTAL_RED_RAMP_START_NM} "
            f"min:{PEDESTAL_FRACTION_AT_MIN_WAVELENGTH} "
            f"max:{PEDESTAL_FRACTION_AT_MAX_WAVELENGTH}"
        ),
    ]

    (
        OUTPUT_PATH
        / "properties"
    ).write_text(
        "\n".join(properties)
        + "\n",
        encoding="utf-8",
    )


# ======================================================================
# Main
# ======================================================================

def main() -> None:

    print("Input bands:")

    for (
        band,
        weight,
        rgb,
        pedestal_fraction,
    ) in zip(
        BANDS,
        WEIGHTS,
        RGB_BASIS,
        BAND_PEDESTAL_FRACTIONS,
    ):

        pixel_cut = property_float_pair(
            band["properties"],
            "hips_pixel_cut",
        )

        pixel_cut_text = (
            "missing"
            if pixel_cut is None
            else (
                f"{pixel_cut[0]:g}"
                f".."
                f"{pixel_cut[1]:g}"
            )
        )

        print(
            f'  {band["id"]:>4}  '
            f'{band["wavelength_nm"]:>4} nm  '
            f"weight={weight:g}  "
            f"pedestal={pedestal_fraction:g}  "
            f"cut={pixel_cut_text}  "
            f"rgb=("
            f"{rgb[0]:.3f},"
            f"{rgb[1]:.3f},"
            f"{rgb[2]:.3f}"
            f")"
        )

    moc = build_union_moc()

    end_order = int(
        moc.max_order
        if END_ORDER is None
        else min(
            END_ORDER,
            moc.max_order,
        )
    )
    end_order = 10

    print(
        f"MOC max_order: {moc.max_order}"
    )

    print(
        f"MOC sky_fraction: {moc.sky_fraction}"
    )

    print(
        f"Processing orders "
        f"{START_ORDER}..{end_order}"
    )

    print(
        f"Output path: {OUTPUT_PATH}"
    )

    write_output_metadata(
        moc,
        end_order,
    )

    for order in range(
        START_ORDER,
        end_order + 1,
    ):

        print(
            f"PROCESSING ORDER {order}"
        )

        params = [
            (
                order,
                int(ipix),
            )
            for ipix in (
                moc
                .degrade_to_order(order)
                .flatten()
            )
        ]

        if (
            order == 3
            and has_input_allsky()
        ):
            params.append(
                (None, None)
            )

        with Pool(
            POOL_SIZE
        ) as pool:

            pool.map(
                process,
                params,
            )


if __name__ == "__main__":
    main()
