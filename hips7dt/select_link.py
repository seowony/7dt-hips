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

ZP: Hipsgen cannot rescale individual files (tested: neither a .hhh
header with expTime= nor skyVal=auto does it), so frames with different
zero points are averaged as they are. zp_report.csv shows the spread per
band; if it is more than ~0.05-0.1 mag, use --scaled-copy ZPREF, which
writes float32 copies scaled to ZPREF (data * 10**(0.4*(ZPREF - zp)))
instead of links - costs disk, so check the size it prints first.
--sky-sub also subtracts each frame's SKYVAL before scaling.

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
    if a.zp_clip > 0:
        # Frames whose ZP is far from the band's median: clouds (ZP low) or a
        # failed calibration (e.g. m625 ZP 29.9 vs median 26.2). Median over
        # all still-eligible frames of the band, all tiles and units.
        cand = df["reason"] == ""
        med = df[cand].groupby("band")["zp"].median()
        off = (df["zp"] - df["band"].map(med)).abs()
        drop(cand & off.notna() & (off > a.zp_clip), "zp_outlier")

    ok = df[df["reason"] == ""].copy()
    ok["_s"] = ok["seeing"].fillna(99.0)
    ok["_d"] = -ok["depth"].fillna(-99.0)
    ok = ok.sort_values(["tile", "band", "_s", "_d", "path"])
    ok["rank"] = ok.groupby(["tile", "band"]).cumcount() + 1
    df["rank"] = ok["rank"]
    # exposure available per tile x band after the quality cuts (before best-N)
    df["avail_exptime"] = ok.groupby(["tile", "band"])["exptime"].transform("sum")
    target = getattr(a, "target_exptime", 0) or 0
    if target > 0:
        # Best-seeing frames until the summed exposure reaches the target, so
        # every band (and tile) gets about the same total exposure. The frame
        # that crosses the target is kept.
        before = ok.groupby(["tile", "band"])["exptime"].cumsum() - ok["exptime"]
        df.loc[ok.index[before >= target], "reason"] = "target_exptime_reached"
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


def exptime_report(df: pd.DataFrame) -> pd.DataFrame:
    """Per tile x band: frames and summed exposure available after the
    quality cuts, and selected."""
    cand = df[df["avail_exptime"].notna()]
    g = cand.groupby(["tile", "band"])
    out = pd.DataFrame({"n_avail": g.size(), "avail_exptime": g["avail_exptime"].first()})
    sel = df[df["selected"] == 1].groupby(["tile", "band"])
    out["n_sel"] = sel.size()
    out["sel_exptime"] = sel["exptime"].sum()
    return out.fillna(0).reset_index()


def zp_by_unit(sel: pd.DataFrame) -> pd.DataFrame:
    g = sel.groupby(["band", "unit"])["zp"]
    out = pd.DataFrame({"n": g.size(), "zp_median": g.median(), "zp_std": g.std()})
    out["offset_vs_band"] = out["zp_median"] - sel.groupby("band")["zp"].median().reindex(
        out.index.get_level_values("band")).to_numpy()
    return out.reset_index()


def write_scaled(src: str, dst: Path, zp: float, zpref: float, sky: float | None):
    """float32 copy of one frame: (data - sky) * 10**(0.4*(zpref - zp)).

    Hipsgen cannot scale individual files: tested on m650 (2026-11-25), a
    HIPSSCAL card in a <name>.hhh header with expTime=HIPSSCAL left the
    pixels unchanged (ratio 1.00), and skyVal=auto multiplied everything by
    ~2.25 instead of removing the sky. So the scaling is done here.
    sky = the frame's own SKYVAL (header) when --sky-sub, else None.
    Written to <dst>.part first and renamed, so an interrupted run never
    leaves a truncated file that a rerun would skip.
    """
    with fits.open(src, memmap=True) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        hdr = hdu.header.copy()
        data = hdu.data.astype(np.float32)
    if sky is not None:
        data -= np.float32(sky)
    scale = np.float32(10 ** (0.4 * (zpref - zp)))
    data *= scale
    for k in ("BSCALE", "BZERO", "BLANK", "ZP_AUTO"):
        hdr.pop(k, None)
    hdr["HIPSZP"] = (zpref, "7dt-hips: pixels scaled to this AB zero point")
    hdr["HIPSSCL"] = (float(scale), "7dt-hips: flux scale applied")
    hdr["HIPSSKY"] = (float(sky) if sky is not None else 0.0, "7dt-hips: sky subtracted [ADU]")
    hdr["HISTORY"] = f"7dt-hips scaled copy of {src}"
    tmp = dst.with_name(dst.name + ".part")
    fits.writeto(tmp, data, hdr, overwrite=True)
    os.replace(tmp, dst)


def _scaled_job(args):
    src, dst, zp, zpref, sky = args
    try:
        write_scaled(src, Path(dst), zp, zpref, sky)
        return ""
    except Exception as exc:  # noqa: BLE001
        return f"{src}: {exc}"


def build_tree(sel: pd.DataFrame, out: Path, scaled_zp: float | None, dry: bool,
               sky_sub: bool = False, workers: int = 1):
    jobs, n_new = [], 0
    for r in sel.itertuples():
        dst = out / r.band / link_name(r)
        if dst.name.endswith(".fz"):
            dst = dst.with_name(dst.name[:-3])
        if dry:
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists() or dst.is_symlink():
            continue
        if scaled_zp is None:
            os.symlink(os.path.abspath(r.path), dst)
            n_new += 1
        else:
            sky = float(r.sky) if sky_sub and pd.notna(r.sky) else None
            if sky_sub and sky is None:
                sky = _header_sky(r.path)
            jobs.append((r.path, str(dst), float(r.zp), scaled_zp, sky))
    if jobs:
        from multiprocessing import Pool
        with Pool(max(1, workers)) as pool:
            errs = [e for e in pool.imap_unordered(_scaled_job, jobs) if e]
        for e in errs:
            print("FAILED", e, file=sys.stderr)
        n_new += len(jobs) - len(errs)
    return n_new


def _header_sky(path: str) -> float | None:
    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdul:
        for h in hdul:
            if h.header.get("NAXIS", 0) >= 2 and "SKYVAL" in h.header:
                return float(h.header["SKYVAL"])
    return None


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
    zp_by_unit(sel).to_csv(out / "zp_by_unit.csv", index=False)
    exp = exptime_report(df)
    exp.to_csv(out / "exptime_report.csv", index=False)
    per_band = exp.groupby("band").agg(sel_min=("sel_exptime", "min"),
                                       sel_max=("sel_exptime", "max"),
                                       avail_min=("avail_exptime", "min"))
    print("summed exposure per tile x band [s]: selected min/max, available min")
    print(per_band.round(0).astype(int).to_string())
    lim = exp["avail_exptime"].min()
    print(f"largest exposure every tile x band can reach: {lim:.0f} s "
          f"(use --target-exptime <= this for a matched depth)")

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
        npix = (sel["naxis1"].fillna(9576) * sel["naxis2"].fillna(6388))  # v2 single size
        gb = npix.sum() * 4 / 1e9
        print(f"--scaled-copy writes ~{gb:.1f} GB")
    n = build_tree(sel, out, a.scaled_copy, a.dry_run, a.sky_sub, a.workers)
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
            dry_run=False, zp_clip=0.3, sky_sub=False, workers=1)
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
        # ZP outlier: give one frame a ZP 2 mag off -> dropped
        idx_df = pd.read_csv(idx)
        idx_df.loc[0, "zp"] = idx_df["zp"].median() + 2.0
        idx2 = Path(tmp) / "idx2.csv"
        idx_df.to_csv(idx2, index=False)
        a2 = argparse.Namespace(**{**vars(a), "index": str(idx2), "out": str(Path(tmp) / "t2")})
        run(a2)
        m2 = pd.read_csv(Path(tmp) / "t2/manifest.csv", dtype=str)
        checks["ZP outlier dropped"] = (m2["reason"] == "zp_outlier").sum() == 1
        # exposure-matched selection: 3 nights x 2 units x 100 s per tile-band
        a3 = argparse.Namespace(**{**vars(a), "out": str(Path(tmp) / "t3"),
                                   "per_tile_band": 0, "target_exptime": 250.0,
                                   "exclude_images": []})
        run(a3)
        e3 = pd.read_csv(Path(tmp) / "t3/exptime_report.csv")
        checks["--target-exptime 250: 3 x 100 s per tile-band"] = (
            (e3["n_sel"] == 3).all() and (e3["sel_exptime"] == 300).all()
            and (e3["avail_exptime"] == 600).all())
        sc = Path(tmp) / "scaled"
        sel = man[man["selected"] == "1"].head(1).assign(zp=lambda d: d["zp"].astype(float))
        build_tree(sel, sc, 25.0, False)
        r = sel.iloc[0]
        f = next(sc.glob("*/*.fits"))
        k = 10 ** (0.4 * (25.0 - r["zp"]))
        ratio = np.median(fits.getdata(f)) / np.median(fits.getdata(r["path"]))
        checks["scaled copy factor"] = abs(ratio - k) < 1e-4
        sc2 = Path(tmp) / "scaled_sky"
        build_tree(sel.assign(sky=np.nan), sc2, 25.0, False, sky_sub=True, workers=2)
        f2 = next(sc2.glob("*/*.fits"))
        sky = fits.getheader(r["path"])["SKYVAL"]
        expect = (fits.getdata(r["path"]).astype(np.float32) - sky) * k
        h2 = fits.getheader(f2)
        checks["sky-sub copy = (data - SKYVAL) * scale, header kept"] = (
            np.allclose(fits.getdata(f2), expect, rtol=1e-5, atol=1e-3)
            and h2["HIPSZP"] == 25.0 and "CTYPE1" in h2
            and not list(sc2.glob("*/*.part")))
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
    ap.add_argument("--per-tile-band", type=int, default=10,
                    help="max frames per tile x band (0 = no limit)")
    ap.add_argument("--target-exptime", type=float, default=0,
                    help="summed exposure per tile x band [s]: best-seeing frames "
                         "until this is reached (matched depth); 0 = off")
    ap.add_argument("--zp-warn", type=float, default=0.1)
    ap.add_argument("--zp-clip", type=float, default=0.3,
                    help="drop frames with |zp - band median| > this [mag]; 0 = off")
    ap.add_argument("--scaled-copy", type=float, metavar="ZPREF")
    ap.add_argument("--sky-sub", action="store_true",
                    help="with --scaled-copy: subtract each frame's SKYVAL first")
    ap.add_argument("--workers", type=int, default=4, help="parallel scaled copies")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        return demo()
    if not a.index:
        ap.error("--index is required")
    sys.path.insert(0, str(Path(__file__).parent))
    from paths import resolve
    a.index, a.out = resolve(a.index), resolve(a.out)
    a.bands = [b for b in a.bands.replace(",", " ").split() if b]
    sys.path.insert(0, str(Path(__file__).parent))
    from frame_index import read_list
    a.tiles = read_list(a.tiles) or []
    a.exclude_images = [p for p in a.exclude_images.split(",") if p]
    return run(a)


if __name__ == "__main__":
    sys.exit(main())
