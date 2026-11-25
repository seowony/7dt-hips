#!/usr/bin/env python3
"""Step 4: colour PNG HiPS from the per-band FITS HiPS.

Thomas Boch's (CDS) make_color_hips.py, cleaned up WITHOUT changing the
colour recipe: every default below is the value his script ended up with
(the last of each repeated assignment). The original is kept verbatim in
docs/thomas/make_color_hips_thomas_v1.py; `--demo` checks that this
version writes pixel-identical tiles.

What changed (plumbing only):
  * command line: --in-root/--out/--orders/--workers/--overwrite
    (his ./hips, END_ORDER=10 hard-coded, 16 processes, always overwrite);
  * resume by default: existing PNG tiles are skipped;
  * hips_order in the properties = the last order processed (was "10");
  * PNG rows are always flipped (his PIL fallback did not flip -> upside-
    down tiles whenever imageio was missing);
  * a band missing at the Allsky gets the Allsky's shape (was 512x512);
    tiles where every band is missing are skipped, not written black;
  * failed tiles are logged and counted instead of killing the pool.

    python3 hips7dt/make_color_hips.py --in-root hips/dp0 --out hips/dp0/color \
        --workers 2
    python3 hips7dt/make_color_hips.py --demo
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
import traceback
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.io import fits
from mocpy import MOC
from PIL import Image

# ----------------------------------------------------------------------
# Thomas's recipe (defaults = final values in his script)
# ----------------------------------------------------------------------
TILE_WIDTH = 512
RED_WEIGHT_ROLLOFF_START_NM = 650
RED_WEIGHT_AT_MAX_WAVELENGTH = 1.00
PEDESTAL_FRACTION_AT_MIN_WAVELENGTH = 0.04
PEDESTAL_RED_RAMP_START_NM = 625
PEDESTAL_FRACTION_AT_MAX_WAVELENGTH = 0.04
RGB_CHANNEL_GAINS = np.array([1.09, 0.90, 1.04], dtype=np.float32) / np.float32(1.5)
RGB_BLACKPOINT = np.array([0.032, 0.010, 0.011], dtype=np.float32)
ORDER_BLACKPOINT_START = 8
ORDER_BLACKPOINT_FULL = 10
ORDER_BLACKPOINT_AT_FULL = np.array([0.004, 0.004, 0.004], dtype=np.float32)
PER_BAND_NORM = "hips_pixel_cut_pedestal"
I_MODE = "luma"
Q = 15.0
I0_MODE = 4.5
GAMMA = 1.25
FINAL_DENOM = 12
SATURATION_SHADOW = 0.14
SATURATION_HIGHLIGHT = 0.80
SATURATION_LUMA_LOW = 0.08
SATURATION_LUMA_HIGH = 0.45
LUMA = (0.2126, 0.7152, 0.0722)

# Filled by configure() in the parent before the pool forks.
CFG: dict = {}


# ----------------------------------------------------------------------
# Bands / properties
# ----------------------------------------------------------------------
def read_properties(path: Path) -> dict[str, str]:
    props = {}
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                props[k.strip()] = v.strip()
    return props


def float_pair(props, key):
    try:
        a, b = props[key].split()[:2]
        return float(a), float(b)
    except (KeyError, ValueError):
        return None


def discover_bands(root: Path, only: list[str] | None = None) -> list[dict]:
    bands = []
    for path in root.iterdir():
        m = re.match(r"^m(\d+)$", path.name)
        if not (path.is_dir() and m and (path / "Moc.fits").exists()):
            continue
        if only and path.name not in only:
            continue
        props = read_properties(path / "properties")
        bands.append({"id": path.name, "path": path, "moc_path": path / "Moc.fits",
                      "properties": props, "wavelength_nm": int(m.group(1)),
                      "pixel_cut": float_pair(props, "hips_pixel_cut")})
    if not bands:
        raise RuntimeError(f"No input HiPS (m<nnn>/Moc.fits) in {root}")
    return sorted(bands, key=lambda b: b["wavelength_nm"])


def wavelength_ramp(x, start, end, v0, v1):
    if end <= start:
        return np.full_like(x, v1, dtype=np.float32)
    t = np.clip((x - start) / (end - start), 0.0, 1.0)
    return (v0 + t * (v1 - v0)).astype(np.float32)


def triangular(x, center, half_width):
    return np.clip(1.0 - np.abs(x - center) / half_width, 0.0, 1.0)


def make_rgb_basis(wl):
    basis = np.stack([triangular(wl, 760.0, 190.0),
                      triangular(wl, 575.0, 155.0),
                      triangular(wl, 450.0, 110.0)], axis=1).astype(np.float32)
    basis /= np.maximum(np.max(basis, axis=1, keepdims=True), 1e-6)
    channel_sum = np.sum(basis, axis=0, keepdims=True)
    basis *= np.mean(channel_sum) / np.maximum(channel_sum, 1e-6)
    return basis


def configure(in_root: Path, out: Path, overwrite: bool, only=None):
    bands = discover_bands(in_root, only)
    wl = np.array([b["wavelength_nm"] for b in bands], dtype=np.float32)
    missing = [b["id"] for b in bands if b["pixel_cut"] is None]
    if missing:
        raise ValueError(f"hips_pixel_cut missing in properties of {missing}")
    CFG.update(
        bands=bands, out=out, overwrite=overwrite, wavelengths=wl,
        weights=wavelength_ramp(wl, RED_WEIGHT_ROLLOFF_START_NM, float(wl.max()),
                                1.0, RED_WEIGHT_AT_MAX_WAVELENGTH),
        pedestal=wavelength_ramp(wl, PEDESTAL_RED_RAMP_START_NM, float(wl.max()),
                                 PEDESTAL_FRACTION_AT_MIN_WAVELENGTH,
                                 PEDESTAL_FRACTION_AT_MAX_WAVELENGTH),
        basis=make_rgb_basis(wl))
    return CFG


# ----------------------------------------------------------------------
# FITS I/O
# ----------------------------------------------------------------------
def tile_suffix(order, ipix):
    if order is None:
        return "Norder3/Allsky.fits"
    return f"Norder{order}/Dir{10000 * (ipix // 10000)}/Npix{ipix}.fits"


def load_fits_any(path: Path) -> np.ndarray:
    with fits.open(path, memmap=True) as hdul:
        for key in ("SCI", "IMAGE", "PRIMARY"):
            try:
                if key in hdul and hdul[key].data is not None:
                    return np.asarray(hdul[key].data, dtype=np.float32)
            except Exception:  # noqa: BLE001
                pass
        for hdu in hdul:
            if getattr(hdu, "data", None) is not None and np.asarray(hdu.data).ndim >= 2:
                return np.asarray(hdu.data).astype(np.float32)
    raise ValueError(f"No 2D image in {path}")


# ----------------------------------------------------------------------
# Colour recipe (unchanged arithmetic)
# ----------------------------------------------------------------------
def normalize_band(band, pixel_cut, pedestal_fraction):
    valid = np.isfinite(band)
    band = np.nan_to_num(band.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    high = pixel_cut[1]
    pedestal = pedestal_fraction * high
    denom = high - pedestal
    if not np.isfinite(denom) or denom <= 0:
        normalized = band
    else:
        normalized = np.clip(band - pedestal, 0.0, denom) / denom
    return np.where(valid, normalized, 0.0).astype(np.float32)


def composite(bands):
    weighted = [w * normalize_band(b, band["pixel_cut"], p)
                for w, b, band, p in zip(CFG["weights"], bands, CFG["bands"],
                                         CFG["pedestal"])]
    rgb_lin = np.stack(weighted, axis=-1) @ CFG["basis"]
    rgb_lin = np.clip(np.nan_to_num(rgb_lin, nan=0.0, posinf=0.0, neginf=0.0), 0.0, None)
    rw, gw, bw = LUMA
    i = np.clip(rw * rgb_lin[..., 0] + gw * rgb_lin[..., 1] + bw * rgb_lin[..., 2],
                0.0, None)
    eps = 1e-6
    x = i / (I0_MODE + eps)
    scale = (np.arcsinh(Q * x) / np.maximum(x, eps)).astype(np.float32)
    rgb = rgb_lin * scale[..., None]
    return np.clip(np.nan_to_num(rgb, nan=0.0, posinf=0.0, neginf=0.0), 0.0, None)


def order_blackpoint(order):
    if order is None or ORDER_BLACKPOINT_FULL <= ORDER_BLACKPOINT_START:
        return np.zeros(3, dtype=np.float32)
    t = np.clip((order - ORDER_BLACKPOINT_START)
                / (ORDER_BLACKPOINT_FULL - ORDER_BLACKPOINT_START), 0.0, 1.0)
    return (t * ORDER_BLACKPOINT_AT_FULL).astype(np.float32)


def display(rgb, order):
    rgb = np.clip(rgb * RGB_CHANNEL_GAINS - RGB_BLACKPOINT - order_blackpoint(order),
                  0.0, None)
    rgb = np.clip(np.nan_to_num(rgb, nan=0.0, posinf=0.0, neginf=0.0), 0.0, None)
    rgb = rgb / (FINAL_DENOM + 1e-6)
    rgb = np.power(rgb, 1.0 / GAMMA, dtype=np.float32)
    rgb = np.clip(rgb, 0.0, 1.0).astype(np.float32)
    # luminance-dependent saturation
    luma = LUMA[0] * rgb[..., 0] + LUMA[1] * rgb[..., 1] + LUMA[2] * rgb[..., 2]
    t = np.clip((luma - SATURATION_LUMA_LOW)
                / (SATURATION_LUMA_HIGH - SATURATION_LUMA_LOW), 0.0, 1.0)
    t = t * t * (3.0 - 2.0 * t)
    sat = SATURATION_SHADOW + t * (SATURATION_HIGHLIGHT - SATURATION_SHADOW)
    out = luma[..., None] + sat[..., None] * (rgb - luma[..., None])
    return np.clip(out, 0.0, 1.0).astype(np.float32)


def save_png(path: Path, rgb):
    u8 = (np.clip(rgb, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    Image.fromarray(np.ascontiguousarray(np.flipud(u8)), mode="RGB").save(path)


# ----------------------------------------------------------------------
# Tiles
# ----------------------------------------------------------------------
def process(param):
    order, ipix = param
    suffix = tile_suffix(order, ipix)
    out_path = CFG["out"] / suffix.replace(".fits", ".png")
    if out_path.exists() and not CFG["overwrite"]:
        return "skipped"
    try:
        loaded = [load_fits_any(b["path"] / suffix) if (b["path"] / suffix).exists()
                  else None for b in CFG["bands"]]
        present = [a for a in loaded if a is not None]
        if not present:
            return "empty"
        shape = present[0].shape
        bands = [a if a is not None else np.full(shape, np.nan, np.float32)
                 for a in loaded]
        rgb = display(composite(bands), order)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        save_png(out_path, rgb)
        return "written"
    except Exception:  # noqa: BLE001
        return "FAILED " + suffix + "\n" + traceback.format_exc(limit=3)


# ----------------------------------------------------------------------
# MOC / metadata
# ----------------------------------------------------------------------
def build_union_moc():
    moc = None
    for b in CFG["bands"]:
        m = MOC.from_fits(b["moc_path"])
        moc = m if moc is None else moc.union(m)
    return moc


def input_max_order():
    orders = [int(b["properties"]["hips_order"]) for b in CFG["bands"]
              if b["properties"].get("hips_order", "").strip().isdigit()]
    return max(orders) if orders else None


def write_metadata(moc, start, end, creator_did, title):
    out = CFG["out"]
    out.mkdir(parents=True, exist_ok=True)
    moc.save(str(out / "Moc.fits"), format="fits", overwrite=True)
    wl = ", ".join(f'{b["id"]}={b["wavelength_nm"]}nm' for b in CFG["bands"])
    g = " ".join(f"{x:g}" for x in RGB_CHANNEL_GAINS)
    bp = " ".join(f"{x:g}" for x in RGB_BLACKPOINT)
    obp = " ".join(f"{x:g}" for x in ORDER_BLACKPOINT_AT_FULL)
    lines = [
        f"creator_did          = {creator_did}",
        f"obs_title            = {title}",
        f"obs_description      = RGB composite from {len(CFG['bands'])} 7DT HiPS bands: {wl}",
        f"moc_sky_fraction     = {moc.sky_fraction:.6E}",
        "hips_builder         = 7dt-hips make_color_hips.py (after T. Boch, CDS)",
        "hips_version         = 1.4",
        "hips_frame           = equatorial",
        f"hips_order           = {end}",
        f"hips_order_min       = {start}",
        f"hips_tile_width      = {TILE_WIDTH}",
        "hips_status          = private master clonableOnce",
        "hips_tile_format     = png",
        "dataproduct_type     = image",
        f"color_per_band_norm  = {PER_BAND_NORM}",
        f"color_i_mode         = {I_MODE}",
        f"color_q              = {Q}",
        f"color_i0             = {I0_MODE}",
        f"color_gamma          = {GAMMA}",
        f"color_final_denom    = {FINAL_DENOM}",
        f"color_saturation     = shadow:{SATURATION_SHADOW} highlight:{SATURATION_HIGHLIGHT} "
        f"luma:{SATURATION_LUMA_LOW}..{SATURATION_LUMA_HIGH}",
        f"color_rgb_gains      = {g}",
        f"color_rgb_blackpoint = {bp}",
        f"color_order_blackpoint = start:{ORDER_BLACKPOINT_START} "
        f"full:{ORDER_BLACKPOINT_FULL} {obp}",
        f"color_red_weight     = start_nm:{RED_WEIGHT_ROLLOFF_START_NM} "
        f"max_weight:{RED_WEIGHT_AT_MAX_WAVELENGTH}",
        f"color_pedestal_frac  = start_nm:{PEDESTAL_RED_RAMP_START_NM} "
        f"min:{PEDESTAL_FRACTION_AT_MIN_WAVELENGTH} max:{PEDESTAL_FRACTION_AT_MAX_WAVELENGTH}",
    ]
    (out / "properties").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(a) -> int:
    configure(Path(a.in_root), Path(a.out), a.overwrite,
              [b for b in a.bands.split(",") if b] or None)
    for b, w, rgb, p in zip(CFG["bands"], CFG["weights"], CFG["basis"], CFG["pedestal"]):
        lo, hi = b["pixel_cut"]
        print(f'  {b["id"]:>4} {b["wavelength_nm"]:>4} nm weight={w:g} pedestal={p:g} '
              f'cut={lo:g}..{hi:g} rgb=({rgb[0]:.3f},{rgb[1]:.3f},{rgb[2]:.3f})')
    moc = build_union_moc()
    start, end = a.order_min, a.order_max
    if end is None:
        end = input_max_order() or int(moc.max_order)
    print(f"MOC max_order {moc.max_order}, sky_fraction {moc.sky_fraction:.3e}; "
          f"orders {start}..{end} -> {CFG['out']}")
    write_metadata(moc, start, end, a.creator_did, a.title)
    has_allsky = any((b["path"] / "Norder3/Allsky.fits").exists() for b in CFG["bands"])
    n_fail = 0
    for order in range(start, end + 1):
        params = [(order, int(i)) for i in moc.degrade_to_order(order).flatten()]
        if order == 3 and has_allsky:
            params.append((None, None))
        with Pool(a.workers) as pool:
            res = pool.map(process, params, chunksize=8)
        counts = {}
        for r in res:
            key = r.split()[0]
            counts[key] = counts.get(key, 0) + 1
            if key == "FAILED":
                print(r, file=sys.stderr)
        n_fail += counts.get("FAILED", 0)
        print(f"order {order}: {len(params)} tiles {counts}", flush=True)
    print(f"done, failed tiles: {n_fail}")
    return 1 if n_fail else 0


# ----------------------------------------------------------------------
# Demo: fake 20-band HiPS, regression against Thomas's original script
# ----------------------------------------------------------------------
def make_fake_hips(root: Path, order=3, cells=(100, 101), missing=("m875", 101), seed=3):
    rng = np.random.default_rng(seed)
    w = 64
    yy, xx = np.mgrid[:TILE_WIDTH, :TILE_WIDTH]
    for k, nm in enumerate(range(400, 900, 25)):
        d = root / f"m{nm}"
        d.mkdir(parents=True)
        hi = 50.0 + 5 * k
        (d / "properties").write_text(f"hips_pixel_cut = -1 {hi}\nhips_order = {order}\n")
        MOC.from_healpix_cells(np.array(cells), np.full(len(cells), order),
                               max_depth=order).save(str(d / "Moc.fits"), format="fits")
        for ip in cells:
            if (f"m{nm}", ip) == missing:
                continue
            img = rng.normal(3.0, 1.0, (TILE_WIDTH, TILE_WIDTH)).astype(np.float32)
            img += 400.0 * np.exp(-((xx - 150) ** 2 + (yy - 300) ** 2) / 30.0) * (nm / 600)
            img += 300.0 * np.exp(-((xx - 400) ** 2 + (yy - 100) ** 2) / 30.0) * (600 / nm) ** 3
            img[:5, :5] = np.nan
            p = d / tile_suffix(order, ip)
            p.parent.mkdir(parents=True, exist_ok=True)
            fits.writeto(p, img)
        allsky = rng.normal(3.0, 1.0, (w * 3, w * 3)).astype(np.float32)
        (d / "Norder3").mkdir(exist_ok=True)
        fits.writeto(d / "Norder3/Allsky.fits", allsky)


def demo() -> int:
    import importlib.util
    import os

    orig_file = Path(__file__).resolve().parents[1] / "docs/thomas/make_color_hips_thomas_v1.py"
    checks = {}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        make_fake_hips(tmp / "hips")
        cwd = os.getcwd()
        os.chdir(tmp)  # Thomas's script reads ./hips at import time
        try:
            spec = importlib.util.spec_from_file_location("thomas", orig_file)
            thomas = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(thomas)
            for p in [(3, 100), (3, 101)]:
                thomas.process(p)
        finally:
            os.chdir(cwd)

        a = argparse.Namespace(in_root=str(tmp / "hips"), out=str(tmp / "mine"),
                               overwrite=False, bands="", order_min=0, order_max=None,
                               workers=2, creator_did="ivo://test/P/7DT/color",
                               title="demo")
        rc = run(a)
        checks["run ok"] = rc == 0
        for ip in (100, 101):
            sfx = tile_suffix(3, ip).replace(".fits", ".png")
            t = np.asarray(Image.open(tmp / "hips/color" / sfx))
            m = np.asarray(Image.open(tmp / "mine" / sfx))
            checks[f"Npix{ip} identical to Thomas's (missing band in 101)"] = (
                t.shape == m.shape and int(np.abs(t.astype(int) - m).max()) == 0)
        img = np.asarray(Image.open(tmp / "mine" / tile_suffix(3, 100).replace(".fits", ".png")))
        checks["tile not blank, has colour"] = img.max() > 100 and np.ptp(
            img.reshape(-1, 3).mean(0)) > 0
        props = read_properties(tmp / "mine/properties")
        checks["hips_order = input order (3)"] = props.get("hips_order") == "3"
        checks["allsky written with its own shape"] = Image.open(
            tmp / "mine/Norder3/Allsky.png").size == (192, 192)
        before = (tmp / "mine" / tile_suffix(3, 100).replace(".fits", ".png")).stat().st_mtime_ns
        run(a)
        checks["resume skips existing tiles"] = (
            tmp / "mine" / tile_suffix(3, 100).replace(".fits", ".png")).stat().st_mtime_ns == before
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-root", default="hips")
    ap.add_argument("--out", default="hips/color")
    ap.add_argument("--bands", default="", help="comma subset, default all m<nnn>")
    ap.add_argument("--order-min", type=int, default=0)
    ap.add_argument("--order-max", type=int, default=None,
                    help="default: max hips_order of the inputs")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--creator-did", default="ivo://SNU/P/7DT/RIS-DP0/color")
    ap.add_argument("--title", default="7DT RIS DP0 colour composite (test)")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    return demo() if a.demo else run(a)


if __name__ == "__main__":
    sys.exit(main())
