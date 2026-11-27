#!/usr/bin/env python3
"""Emission-line map HiPS from three per-band FITS HiPS (default: H-alpha).

The broad RGB composite hides what is unique to 7DT: 25-nm bands. This
shows a line band against its two neighbours:

    excess = line - (s1 * cont1 + s2 * cont2) / 2

with every band's sky removed. The scales s_b put the continuum bands on
the line band's flux scale; they are measured from the data as the median
of line/cont over bright (star) pixels of sample tiles, so an average star
gives excess = 0 and no zero points are needed. H-alpha at z ~ 0 (656.3 nm)
falls in m650; default continuum = m625 and m675.

Picture: the continuum as a grey image; pixels with excess above
--thr sigma turn RED (emission), below -thr sigma BLUE (line band fainter:
absorption, or an artefact in a continuum band), saturating at --sat sigma.
sigma = robust noise of the excess on empty sky.

Read it with care:
  * the bands have different seeing, so bright stars leave a red core with
    a blue ring (or the reverse) - a PSF residual, not a line;
  * satellite trails, cosmic rays and ghosts in ONE band show up at full
    strength - red if in the line band, blue if in a continuum band. That
    also makes this a per-band artefact finder;
  * very red or blue stars deviate from the "average star" scaling.
A real emitter is red, has no blue ring, and follows the object's shape.

    python3 hips7dt/make_line_hips.py --in-root hips/dp0_n10 \
        --out hips/dp0_n10/halpha --workers 8
    python3 hips7dt/make_line_hips.py --line m500 --cont m475,m525 ...   # [OIII]
    python3 hips7dt/make_line_hips.py --demo
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import traceback
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.io import fits
from mocpy import MOC
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
import make_color_hips as mc  # noqa: E402

CFG: dict = {}


def sample_tiles(band: dict, n_tiles: int) -> list[Path]:
    order = band["properties"].get("hips_order", "").strip()
    tiles = sorted(band["path"].glob(f"Norder{order}/Dir*/Npix*.fits")) if order else []
    return tiles[:: max(1, len(tiles) // n_tiles)][:n_tiles]


def calibrate(line: dict, conts: list[dict], n_tiles: int = 40, bright: float = 30.0):
    """Sky + noise per band, continuum scales from bright pixels, noise of
    the excess and of the continuum."""
    sky, sig = {}, {}
    for b in [line] + conts:
        sky[b["id"]], sig[b["id"]] = mc.measure_sky(b, n_tiles)
    ratios = {c["id"]: [] for c in conts}
    tiles = sample_tiles(line, n_tiles)
    loaded = []
    for t in tiles:
        rel = t.relative_to(line["path"])
        if not all((c["path"] / rel).exists() for c in conts):
            continue
        ln = mc.load_fits_any(t) - sky[line["id"]]
        cs = [mc.load_fits_any(c["path"] / rel) - sky[c["id"]] for c in conts]
        loaded.append((ln, cs))
        for c, cd in zip(conts, cs):
            m = (np.isfinite(ln) & np.isfinite(cd) & (ln > bright * sig[line["id"]])
                 & (cd > bright * sig[c["id"]]))
            if m.sum() > 20:
                ratios[c["id"]].append(float(np.median(ln[m] / cd[m])))
    if not loaded:
        raise RuntimeError("no tile is present in all three bands")
    scale = {}
    for c in conts:
        if not ratios[c["id"]]:
            raise RuntimeError(f"no bright pixels to scale {c['id']} to {line['id']}")
        scale[c["id"]] = float(np.median(ratios[c["id"]]))
    exc, cont = [], []
    for ln, cs in loaded:
        cm = sum(scale[c["id"]] * cd for c, cd in zip(conts, cs)) / len(conts)
        e = ln - cm
        ok = np.isfinite(e)
        exc.append(e[ok][::17])
        cont.append(cm[ok][::17])
    exc, cont = np.concatenate(exc), np.concatenate(cont)
    mad = lambda x: float(1.4826 * np.median(np.abs(x - np.median(x))))  # noqa: E731
    return dict(sky=sky, sig=sig, scale=scale, sig_exc=mad(exc), sig_cont=mad(cont),
                exc_offset=float(np.median(exc)))


def render(line, conts, cal, thr, sat):
    """line, conts: sky-subtracted arrays (NaN = no data). Returns RGB 0..1."""
    cm = sum(s * c for s, c in zip(cal["scales"], conts)) / len(conts)
    valid = np.isfinite(line) & np.isfinite(cm)
    cm = np.where(valid, cm, 0.0)
    exc = np.where(valid, line - cm - cal["exc_offset"], 0.0)
    grey = np.clip(np.arcsinh(np.clip(cm, 0, None) / (5.0 * cal["sig_cont"]))
                   / np.arcsinh(200.0), 0.0, 1.0) * 0.85
    e = exc / cal["sig_exc"]
    t = np.clip((np.abs(e) - thr) / (sat - thr), 0.0, 1.0)
    pos, neg = t * (e > 0), t * (e < 0)
    r = np.maximum(grey * (1 - neg), pos)
    g = grey * (1 - pos) * (1 - neg)
    b = np.maximum(grey * (1 - pos), neg)
    rgb = np.stack([r, g, b], axis=-1)
    return np.where(valid[..., None], rgb, 0.0).astype(np.float32)


def process(param):
    order, ipix = param
    suffix = mc.tile_suffix(order, ipix)
    out_path = CFG["out"] / suffix.replace(".fits", ".png")
    if out_path.exists() and not CFG["overwrite"]:
        return "skipped"
    try:
        paths = [b["path"] / suffix for b in CFG["bands"]]
        if not all(p.exists() for p in paths):
            return "empty"
        arrs = [mc.load_fits_any(p) - CFG["cal"]["sky"][b["id"]]
                for p, b in zip(paths, CFG["bands"])]
        rgb = render(arrs[0], arrs[1:], CFG["cal"], CFG["thr"], CFG["sat"])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        mc.save_png(out_path, rgb)
        return "written"
    except Exception:  # noqa: BLE001
        return "FAILED " + suffix + "\n" + traceback.format_exc(limit=3)


def run(a) -> int:
    cont_ids = [c for c in a.cont.split(",") if c]
    bands = {b["id"]: b for b in mc.discover_bands(Path(a.in_root), [a.line] + cont_ids)}
    missing = [x for x in [a.line] + cont_ids if x not in bands]
    if missing:
        raise SystemExit(f"band HiPS not found under {a.in_root}: {missing}")
    line, conts = bands[a.line], [bands[c] for c in cont_ids]
    cal = calibrate(line, conts, a.sample_tiles)
    cal["scales"] = [cal["scale"][c["id"]] for c in conts]
    CFG.update(bands=[line] + conts, out=Path(a.out), overwrite=a.overwrite, cal=cal,
               thr=a.thr, sat=a.sat)
    for b in [line] + conts:
        s = "" if b is line else f"  scale to {line['id']} = {cal['scale'][b['id']]:.3f}"
        print(f"  {b['id']}: sky={cal['sky'][b['id']]:.3g} sigma={cal['sig'][b['id']]:.3g}{s}")
    print(f"  excess: sigma={cal['sig_exc']:.3g} offset={cal['exc_offset']:.3g}; "
          f"colour from {a.thr} to {a.sat} sigma")

    moc = MOC.from_fits(line["moc_path"])
    for c in conts:
        moc = moc.intersection(MOC.from_fits(c["moc_path"]))
    end = a.order_max if a.order_max is not None else int(
        line["properties"].get("hips_order", moc.max_order))
    out = CFG["out"]
    out.mkdir(parents=True, exist_ok=True)
    moc.save(str(out / "Moc.fits"), format="fits", overwrite=True)
    (out / "properties").write_text("\n".join([
        f"creator_did          = {a.creator_did}",
        f"obs_title            = 7DT line map: {a.line} minus {a.cont}",
        "obs_description      = grey = continuum; red = line-band excess, blue = deficit",
        "hips_builder         = 7dt-hips make_line_hips.py",
        "hips_version         = 1.4", "hips_frame           = equatorial",
        f"hips_order           = {end}", f"hips_order_min       = {a.order_min}",
        f"hips_tile_width      = {mc.TILE_WIDTH}",
        "hips_status          = private master clonableOnce",
        "hips_tile_format     = png", "dataproduct_type     = image",
        f"line_band            = {a.line}", f"line_continuum       = {a.cont}",
        "line_cont_scales     = " + " ".join(f"{s:.4f}" for s in cal["scales"]),
        f"line_sigma_excess    = {cal['sig_exc']:.4g}",
        f"line_threshold_sigma = {a.thr} {a.sat}",
    ]) + "\n", encoding="utf-8")

    n_fail = 0
    for order in range(a.order_min, end + 1):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            cells = moc.degrade_to_order(order).flatten()
        with Pool(a.workers) as pool:
            res = pool.map(process, [(order, int(i)) for i in cells], chunksize=8)
        counts = {}
        for r in res:
            k = r.split()[0]
            counts[k] = counts.get(k, 0) + 1
            if k == "FAILED":
                print(r, file=sys.stderr)
        n_fail += counts.get("FAILED", 0)
        print(f"order {order}: {len(cells)} tiles {counts}", flush=True)
    print(f"done, failed tiles: {n_fail} -> {out}")
    return 1 if n_fail else 0


def demo() -> int:
    """Three fake bands on different flux scales; stars, one emitter, one
    absorber, one single-band artefact."""
    rng = np.random.default_rng(5)
    w, order, cell = 256, 3, 100
    yy, xx = np.mgrid[:w, :w]
    blob = lambda x, y, a, s=6.0: a * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / s)  # noqa: E731
    true_scale = {"m625": 1.0 / 0.785, "m650": 1.0, "m675": 1.0 / 1.163}  # counts vs m650
    stars = [(rng.uniform(10, w - 10), rng.uniform(10, w - 10), rng.uniform(100, 2000))
             for _ in range(40)]                       # same stars in every band
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "hips"
        for band, k in true_scale.items():
            d = root / band
            (d / f"Norder{order}/Dir0").mkdir(parents=True)
            img = rng.normal(2.0, 1.0, (w, w))
            for sx, sy, amp in stars:                  # ordinary stars
                img += blob(sx, sy, amp) * k
            img += blob(60, 60, 300) * k               # emitter: continuum ...
            img += blob(200, 200, 300) * k             # absorber: continuum ...
            if band == "m650":
                img += blob(60, 60, 400)               # ... + line
                img -= blob(200, 200, 200)             # ... - line
            if band == "m675":
                img[120:123, 20:120] += 80             # trail in one continuum band
            fits.writeto(d / f"Norder{order}/Dir0/Npix{cell}.fits", img.astype(np.float32))
            (d / "properties").write_text(f"hips_order = {order}\nhips_pixel_cut = 0 50\n")
            MOC.from_healpix_cells(np.array([cell]), np.array([order]),
                                   max_depth=order).save(str(d / "Moc.fits"), format="fits")
        a = argparse.Namespace(in_root=str(root), out=str(Path(tmp) / "ha"), line="m650",
                               cont="m625,m675", thr=3.0, sat=15.0, sample_tiles=40,
                               order_min=3, order_max=3, workers=1, overwrite=False,
                               creator_did="ivo://test/ha")
        rc = run(a)
        im = np.asarray(Image.open(Path(tmp) / "ha" / f"Norder{order}/Dir0/Npix{cell}.png")
                        ).astype(float)[::-1]          # undo the PNG flip
        cal = CFG["cal"]
        r, g, b = im[..., 0], im[..., 1], im[..., 2]
        star = np.unravel_index(np.argmax(g), g.shape)   # brightest grey pixel
        empty = (np.abs(r - g) < 2) & (np.abs(b - g) < 2)
        checks = {
            "run ok": rc == 0,
            "continuum scales recovered (0.785, 1.163)":
                abs(cal["scale"]["m625"] - 0.785) < 0.02 and abs(cal["scale"]["m675"] - 1.163) < 0.03,
            "emitter is red": r[60, 60] > 200 and g[60, 60] < 60 and b[60, 60] < 60,
            "absorber is blue": b[200, 200] > 200 and r[200, 200] < 60,
            "ordinary star stays grey": abs(r[star] - g[star]) < 15 and abs(b[star] - g[star]) < 15,
            "continuum-band trail shows blue": b[121, 60] > 150 and r[121, 60] < 60,
            "> 97% of pixels uncoloured": empty.mean() > 0.97,
        }
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-root", default="hips")
    ap.add_argument("--out", default="hips/halpha")
    ap.add_argument("--line", default="m650", help="band holding the line (H-alpha: m650)")
    ap.add_argument("--cont", default="m625,m675", help="continuum bands")
    ap.add_argument("--thr", type=float, default=3.0, help="colour starts at this many sigma")
    ap.add_argument("--sat", type=float, default=15.0, help="full colour at this many sigma")
    ap.add_argument("--sample-tiles", type=int, default=40)
    ap.add_argument("--order-min", type=int, default=0)
    ap.add_argument("--order-max", type=int, default=None)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--creator-did", default="ivo://SNU/P/7DT/RIS-DP0/halpha")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        return demo()
    from paths import resolve
    a.in_root, a.out = resolve(a.in_root), resolve(a.out)
    return run(a)


if __name__ == "__main__":
    sys.exit(main())
