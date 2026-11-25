#!/usr/bin/env python3
"""Step 2: choose frames per (tile, band) and build the Hipsgen input tree.

Reads the frame index (frame_index.py) and writes

    <out>/<band>/<tile>_<night>_<unit>__<original name>  ->  /lyman/...  (symlink)
    <out>/manifest.csv                                      (what was chosen, why)
    <out>/zp_report.csv                                     (ZP spread per band)

No pixel data is copied. Selection per (tile, band):
  1. frame readable, WCS present, band in --bands
  2. not excluded (--exclude-dates, --exclude-images CSV with
     tile,obs_date,band columns, e.g. the RIS audit lists); a frame is
     excluded when its UT date OR its night matches
  3. seeing <= --max-seeing (when seeing is known)
  4. best --per-tile-band frames by seeing (then depth); 0 = keep all

ZP: Hipsgen cannot rescale individual files, so frames with different
zero points are averaged as they are. zp_report.csv shows the spread per
band; if it is more than ~0.05-0.1 mag, use --scaled-copy ZPREF, which
writes float32 copies scaled to ZPREF (data * 10**(0.4*(ZPREF - zp)))
instead of links - costs disk, so check the size it prints first.

    python3 hips7dt/select_link.py --index index/dp0_frames.csv \
        --out hips_in/dp0 --per-tile-band 10 --dry-run
    python3 hips7dt/select_link.py --demo
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time

MEDIUM_BANDS = [f"m{w}" for w in range(400, 900, 25)]


def load_exclusions(dates: str | None, images_csv: list[str]):
    bad_dates = set((dates or "").replace(",", " ").split())
    bad_images = set()
    for path in images_csv:
        df = pd.read_csv(path, dtype=str)
        for _, r in df.iterrows():
            bad_images.add((r["tile"], str(r["obs_date"])[:10], r["band"]))
    return bad_dates, bad_images


def select(df: pd.DataFrame, a) -> pd.DataFrame:
    df = df.copy()
    df["reason"] = ""
    df["ut_date"] = [Time(m, format="mjd").to_value("iso", subfmt="date")
                     if pd.notna(m) else "" for m in df["mjd"]]
    df["night"] = df["night"].fillna("").astype(str)

    def drop(mask, why):
        df.loc[mask & (df["reason"] == ""), "reason"] = why

    drop(df["error"].fillna("") != "", "unreadable")
    drop(df["has_wcs"].fillna(0).astype(int) != 1, "no_wcs")
    drop(~df["band"].isin(a.bands), "band_not_used")
    if a.tiles:
        drop(~df["tile"].isin(a.tiles), "tile_not_used")
    bad_dates, bad_images = load_exclusions(a.exclude_dates, a.exclude_images)
    drop(df["ut_date"].isin(bad_dates) | df["night"].isin(bad_dates), "excluded_date")
    if bad_images:
        hit = [((t, u, b) in bad_images) or ((t, n, b) in bad_images)
               for t, u, n, b in zip(df["tile"], df["ut_date"], df["night"], df["band"])]
        drop(pd.Series(hit, index=df.index), "excluded_image")
    if a.require_zp:
        drop(df["zp"].isna(), "no_zp")
    if a.max_seeing:
        drop(df["seeing"].notna() & (df["seeing"] > a.max_seeing), "seeing")

    ok = df[df["reason"] == ""].copy()
    ok["_s"] = ok["seeing"].fillna(99.0)
    ok["_d"] = -ok["depth"].fillna(-99.0)
    ok = ok.sort_values(["tile", "band", "_s", "_d", "path"])
    ok["rank"] = ok.groupby(["tile", "band"]).cumcount() + 1
    df["rank"] = ok["rank"]
    if a.per_tile_band > 0:
        drop(df["rank"].fillna(1e9) > a.per_tile_band, "not_in_best_n")
    df["selected"] = (df["reason"] == "").astype(int)
    return df.drop(columns=["_s", "_d"], errors="ignore")


def link_name(r) -> str:
    return f"{r.tile}_{r.night or r.ut_date}_{r.unit}__{os.path.basename(r.path)}"


def zp_report(sel: pd.DataFrame) -> pd.DataFrame:
    g = sel.groupby("band")["zp"]
    rep = pd.DataFrame({
        "n": sel.groupby("band").size(),
        "zp_median": g.median(), "zp_std": g.std(),
        "zp_min": g.min(), "zp_max": g.max(),
        "n_units": sel.groupby("band")["unit"].nunique(),
        "seeing_median": sel.groupby("band")["seeing"].median(),
    })
    rep["zp_range"] = rep["zp_max"] - rep["zp_min"]
    return rep.reset_index()


def build_tree(sel: pd.DataFrame, out: Path, scaled_zp: float | None, dry: bool):
    n_new = 0
    for r in sel.itertuples():
        dst = out / r.band / link_name(r)
        if dry:
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists() or dst.is_symlink():
            continue
        if scaled_zp is None:
            os.symlink(os.path.abspath(r.path), dst)
        else:
            with fits.open(r.path, memmap=True) as hdul:
                hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
                scale = 10 ** (0.4 * (scaled_zp - r.zp))
                hdr = hdu.header.copy()
                for k in ("BSCALE", "BZERO", "BLANK"):
                    hdr.pop(k, None)
                hdr["HISTORY"] = f"scaled x{scale:.6f} to ZP {scaled_zp} from {r.path}"
                fits.writeto(dst.with_suffix(".fits") if str(dst).endswith(".fz") else dst,
                             (hdu.data.astype(np.float32) * scale), hdr)
        n_new += 1
    return n_new


def run(a) -> int:
    df = pd.read_csv(a.index, dtype={"tile": str, "band": str, "unit": str,
                                     "night": str, "error": str})
    df = select(df, a)
    sel = df[df["selected"] == 1]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "manifest.csv", index=False)
    rep = zp_report(sel)
    rep.to_csv(out / "zp_report.csv", index=False)

    print("dropped:", df.loc[df["selected"] == 0, "reason"].value_counts().to_dict())
    print(f"selected {len(sel)} frames, "
          f"{sel.groupby(['tile', 'band']).ngroups} tile-bands")
    with pd.option_context("display.width", 140, "display.max_rows", 50):
        print(rep.round(3).to_string(index=False))
    big = rep[rep["zp_range"] > a.zp_warn]
    if len(big):
        print(f"WARNING: ZP range > {a.zp_warn} mag in {list(big['band'])} - "
              "frames will be averaged with different scales (seams); "
              "consider --scaled-copy.")
    if a.scaled_copy is not None:
        gb = (sel["naxis1"] * sel["naxis2"]).sum() * 4 / 1e9
        print(f"--scaled-copy writes ~{gb:.1f} GB")
    n = build_tree(sel, out, a.scaled_copy, a.dry_run)
    print(("dry run: nothing written to the tree" if a.dry_run
           else f"{n} new entries under {out}/<band>/"))
    return 0


def demo() -> int:
    sys.path.insert(0, str(Path(__file__).parent))
    import frame_index as fi

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "lyman"
        fi.make_demo_frames(root, ["T02524", "T02666"], ["m400", "m650", "r"],
                            ["20251110", "20251111", "20251112"],
                            ["7DT01", "7DT02"])
        rows = [fi.index_one(f) for f in
                fi.find_files([str(root / "**" / "*.fits")], None)]
        idx = Path(tmp) / "idx.csv"
        fi.write_csv(rows, idx)
        bad = Path(tmp) / "bad.csv"
        bad.write_text("release,tile,obs_date,band\ndp0,T02524,2025-11-11,m650\n")
        out = Path(tmp) / "tree"
        a = argparse.Namespace(
            index=str(idx), out=str(out), bands=MEDIUM_BANDS, tiles=None,
            exclude_dates="", exclude_images=[str(bad)], require_zp=True,
            max_seeing=0, per_tile_band=3, zp_warn=0.1, scaled_copy=None,
            dry_run=False)
        run(a)
        man = pd.read_csv(out / "manifest.csv", dtype=str)
        links = sorted((out).glob("*/*.fits"))
        n_t2524_m650 = len(list((out / "m650").glob("T02524_*")))
        checks = {
            "broadband r not used": not (out / "r").exists(),
            "<= 3 per tile-band": all(len(list((out / b).glob(f"{t}_*"))) <= 3
                                      for b in ("m400", "m650")
                                      for t in ("T02524", "T02666")),
            "excluded image dropped": n_t2524_m650 == 3 and
                (man["reason"] == "excluded_image").sum() == 2,
            "entries are symlinks to the originals": links and all(
                p.is_symlink() and Path(os.readlink(p)).exists() for p in links),
            "rerun adds nothing": build_tree(
                man[man["selected"] == "1"].assign(
                    zp=lambda d: d["zp"].astype(float)), out, None, False) == 0,
        }
        sc = Path(tmp) / "scaled"
        sel = man[man["selected"] == "1"].head(1).assign(zp=lambda d: d["zp"].astype(float))
        build_tree(sel, sc, 25.0, False)
        r = sel.iloc[0]
        f = next(sc.glob("*/*.fits"))
        ratio = np.median(fits.getdata(f)) / np.median(fits.getdata(r["path"]))
        checks["scaled copy factor"] = abs(ratio - 10 ** (0.4 * (25.0 - r["zp"]))) < 1e-4
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index")
    ap.add_argument("--out", default="hips_in")
    ap.add_argument("--bands", default=",".join(MEDIUM_BANDS))
    ap.add_argument("--tiles", default="", help="comma list or file; default all")
    ap.add_argument("--exclude-dates", default="", help="YYYY-MM-DD list")
    ap.add_argument("--exclude-images", default="",
                    help="comma list of CSVs with tile,obs_date,band")
    ap.add_argument("--require-zp", action="store_true")
    ap.add_argument("--max-seeing", type=float, default=0.0, help="arcsec, 0 = off")
    ap.add_argument("--per-tile-band", type=int, default=10)
    ap.add_argument("--zp-warn", type=float, default=0.1)
    ap.add_argument("--scaled-copy", type=float, metavar="ZPREF")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        return demo()
    if not a.index:
        ap.error("--index is required")
    a.bands = [b for b in a.bands.replace(",", " ").split() if b]
    if a.tiles and os.path.exists(a.tiles):
        a.tiles = Path(a.tiles).read_text().split()
    else:
        a.tiles = [t for t in a.tiles.split(",") if t]
    a.exclude_images = [p for p in a.exclude_images.split(",") if p]
    return run(a)


if __name__ == "__main__":
    sys.exit(main())
