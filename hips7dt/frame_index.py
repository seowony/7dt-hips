#!/usr/bin/env python3
"""Step 1: header-only index of 7DT single processed science frames.

Reads ONLY headers (no pixel data), so it is cheap even over NFS.
One CSV row per frame: path, tile, band, unit, night, mjd, exptime,
zp, seeing, depth, sky, image centre (from the WCS), has_wcs.

Header keyword names of the v2 frames are not known yet, so every
quantity is looked up from a candidate list (first key present wins)
and the key actually used is recorded (zp_key, seeing_key, ...).
Check them first on one real frame:

    python3 hips7dt/frame_index.py --dump-header /lyman/.../one_frame.fits

then run the index (DP0 tiles):

    python3 hips7dt/frame_index.py \
        --glob '/lyman/data2/<v2 root>/*/{tile}/*/*.fits' \
        --tiles config/dp0_tiles.txt --out index/dp0_frames.csv --workers 2

`{tile}` in --glob is expanded per tile (keeps the scan to those tiles).
Paths containing /coadd/ are skipped (single frames only).

    python3 hips7dt/frame_index.py --demo     # no data needed
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

KEYS = {
    "band": ["FILTER", "FILTNAME"],
    "tile": ["OBJECT", "TILE"],
    "unit": ["TELESCOP", "INSTRUME", "OBSERVAT"],
    "mjd": ["MJD-OBS", "MJD"],
    "date_obs": ["DATE-OBS"],
    "exptime": ["EXPTIME", "EXPOSURE"],
    "zp": ["ZP_AUTO", "ZP", "ZP5_1", "ZP_5", "PHOTZP", "MAGZERO", "MAGZP"],
    "seeing": ["SEEING", "FWHM_ARC", "FWHM", "PEEING"],
    "depth": ["UL5_AUTO", "UL5_5", "UL5_1", "DEPTH", "LIMMAG"],
    "sky": ["SKYVAL", "SKY", "BACKGND"],
    "skysig": ["SKYSIG", "BKGRMS", "SKYRMS"],
}

TILE_RE = re.compile(r"(T\d{5})")
BAND_RE = re.compile(r"(?<![a-z0-9])(m\d{3}|[ugriz])(?![a-z0-9])")
UNIT_RE = re.compile(r"(7DT\d{2})")
FIELDS = (
    ["path", "tile", "band", "unit", "night", "mjd", "date_obs", "exptime",
     "zp", "seeing", "depth", "sky", "skysig", "ra_cen", "dec_cen",
     "naxis1", "naxis2", "has_wcs"]
    + [f"{k}_key" for k in ("zp", "seeing", "depth", "sky")]
    + ["error"]
)


def image_hdu(path: str) -> int:
    """First HDU with a 2D image (0 for plain FITS, 1 for .fz)."""
    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdul:
        for i, hdu in enumerate(hdul):
            if hdu.header.get("NAXIS", 0) >= 2 or hdu.header.get("ZNAXIS", 0) >= 2:
                return i
    return 0


def first(header, kind):
    for key in KEYS[kind]:
        if key in header and header[key] not in (None, ""):
            return header[key], key
    return None, ""


def as_float(v):
    try:
        x = float(v)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def index_one(path: str) -> dict:
    row = {f: "" for f in FIELDS}
    row["path"] = path
    try:
        hdr = fits.getheader(path, ext=image_hdu(path))
        name = os.path.basename(path)

        tile, _ = first(hdr, "tile")
        m = TILE_RE.search(str(tile or "")) or TILE_RE.search(path)
        row["tile"] = m.group(1) if m else ""

        band, _ = first(hdr, "band")
        band = str(band or "").strip()
        if not band:
            m = BAND_RE.search(name)
            band = m.group(1) if m else ""
        row["band"] = band

        unit, _ = first(hdr, "unit")
        m = UNIT_RE.search(str(unit or "")) or UNIT_RE.search(path)
        row["unit"] = m.group(1) if m else str(unit or "")

        for kind in ("mjd", "exptime"):
            row[kind] = as_float(first(hdr, kind)[0])
        row["date_obs"] = str(first(hdr, "date_obs")[0] or "")
        for kind in ("zp", "seeing", "depth", "sky", "skysig"):
            val, key = first(hdr, kind)
            row[kind] = as_float(val)
            if f"{kind}_key" in row:
                row[f"{kind}_key"] = key

        # Night = local date of the evening (UT - 12 h); site-agnostic enough
        # to group the frames of one night together.
        if row["mjd"] is not None:
            from astropy.time import Time
            row["night"] = Time(row["mjd"] - 0.5, format="mjd").to_value(
                "iso", subfmt="date")

        nx = int(hdr.get("ZNAXIS1", hdr.get("NAXIS1", 0)))
        ny = int(hdr.get("ZNAXIS2", hdr.get("NAXIS2", 0)))
        row["naxis1"], row["naxis2"] = nx, ny
        try:
            w = WCS(hdr, naxis=2)
            ok = w.has_celestial and nx > 0 and ny > 0
            if ok:
                ra, dec = w.all_pix2world([[nx / 2 + 0.5, ny / 2 + 0.5]], 1)[0]
                row["ra_cen"], row["dec_cen"] = float(ra), float(dec)
            row["has_wcs"] = int(bool(ok))
        except Exception as exc:  # noqa: BLE001
            row["has_wcs"] = 0
            row["error"] = f"wcs: {exc}"[:200]
    except Exception as exc:  # noqa: BLE001
        row["error"] = str(exc)[:200]
    return row


def find_files(patterns: list[str], tiles: list[str] | None) -> list[str]:
    files = set()
    for pat in patterns:
        expanded = ([pat.replace("{tile}", t) for t in tiles]
                    if tiles and "{tile}" in pat else [pat])
        for p in expanded:
            files.update(glob.glob(p, recursive=True))
    return sorted(f for f in files
                  if "/coadd/" not in f and f.endswith((".fits", ".fits.fz", ".fit")))


def read_list(path_or_csv: str | None) -> list[str] | None:
    if not path_or_csv:
        return None
    if os.path.exists(path_or_csv):
        text = Path(path_or_csv).read_text().replace(",", " ")
        return [w for w in text.split() if not w.startswith("#")]
    return [w for w in path_or_csv.replace(",", " ").split() if w]


def write_csv(rows, out: Path):
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=FIELDS)
        wr.writeheader()
        for r in rows:
            wr.writerow({k: ("" if r[k] is None else r[k]) for k in FIELDS})


def summarise(rows):
    ok = [r for r in rows if not r["error"]]
    print(f"frames: {len(rows)}  readable: {len(ok)}  "
          f"with WCS: {sum(1 for r in ok if r['has_wcs'])}")
    for kind in ("zp", "seeing", "depth", "sky"):
        keys = {}
        for r in ok:
            keys[r[f"{kind}_key"] or "(missing)"] = keys.get(
                r[f"{kind}_key"] or "(missing)", 0) + 1
        print(f"  {kind:7s} keys used: {keys}")
    bands = sorted({r["band"] for r in ok})
    tiles = sorted({r["tile"] for r in ok})
    print(f"  tiles: {tiles}")
    print(f"  bands: {bands}")


def make_demo_frames(root: Path, tiles, bands, nights, units, seed=1):
    """Small fake frames with a TAN WCS and 7DT-like headers."""
    rng = np.random.default_rng(seed)
    paths = []
    for ti, tile in enumerate(tiles):
        ra0, dec0 = 150.0 + 1.3 * ti, -30.0
        for band in bands:
            for ni, night in enumerate(nights):
                for unit in units:
                    d = root / night / tile / band
                    d.mkdir(parents=True, exist_ok=True)
                    hdr = fits.Header()
                    hdr.update({
                        "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
                        "CRVAL1": ra0 + rng.normal(0, 0.01),
                        "CRVAL2": dec0 + rng.normal(0, 0.01),
                        "CRPIX1": 32.5, "CRPIX2": 24.5,
                        "CD1_1": -0.000140, "CD1_2": 0.0,
                        "CD2_1": 0.0, "CD2_2": 0.000140,
                        "FILTER": band, "OBJECT": tile, "TELESCOP": unit,
                        "MJD-OBS": 60990.6 + ni, "EXPTIME": 100.0,
                        "ZP_AUTO": 25.0 + rng.normal(0, 0.05),
                        "SEEING": float(rng.uniform(1.5, 4.0)),
                        "UL5_5": float(rng.uniform(18.5, 20.0)),
                        "SKYVAL": float(rng.uniform(50, 300)),
                    })
                    data = rng.normal(100, 5, (48, 64)).astype(np.float32)
                    p = d / f"calib_{unit}_{tile}_{night}_{band}_100.0s.fits"
                    fits.writeto(p, data, hdr, overwrite=True)
                    paths.append(str(p))
            cd = root / nights[0] / tile / band / "coadd"
            cd.mkdir(exist_ok=True)
            fits.writeto(cd / "skip_me_coadd.fits",
                         np.zeros((4, 4), np.float32), overwrite=True)
    return paths


def demo() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "lyman"
        made = make_demo_frames(root, ["T02524", "T02666"], ["m400", "m650"],
                                ["20251110", "20251111"], ["7DT01", "7DT02"])
        files = find_files([str(root / "*" / "{tile}" / "*" / "*.fits")],
                           ["T02524", "T02666"])
        rows = [index_one(f) for f in files]
        summarise(rows)
        checks = {
            "all frames found, coadd skipped": len(files) == len(made) == 16,
            "no read errors": all(not r["error"] for r in rows),
            "tile/band/unit parsed": all(r["tile"] and r["band"] and
                                         r["unit"].startswith("7DT") for r in rows),
            "zp/seeing via candidate keys": all(r["zp_key"] == "ZP_AUTO" and
                                                r["seeing_key"] == "SEEING" for r in rows),
            "wcs centre": all(r["has_wcs"] == 1 and abs(r["dec_cen"] + 30) < 0.1
                              for r in rows),
        }
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--glob", action="append", default=[],
                    help="file glob; {tile} expands per tile; repeatable")
    ap.add_argument("--tiles", help="file or comma list of tiles")
    ap.add_argument("--out", default="index/frames.csv")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--dump-header", metavar="FITS")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()

    if a.demo:
        return demo()
    if a.dump_header:
        ext = image_hdu(a.dump_header)
        print(f"# HDU {ext}")
        print(repr(fits.getheader(a.dump_header, ext=ext)))
        return 0
    if not a.glob:
        ap.error("--glob is required")

    tiles = read_list(a.tiles)
    files = find_files(a.glob, tiles)
    print(f"{len(files)} files matched")
    if not files:
        return 1
    with Pool(a.workers) as pool:
        rows = pool.map(index_one, files, chunksize=32)
    if tiles:
        rows = [r for r in rows if r["tile"] in tiles or r["error"]]
    write_csv(rows, Path(a.out))
    summarise(rows)
    print(f"written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
