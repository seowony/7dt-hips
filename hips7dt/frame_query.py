#!/usr/bin/env python3
"""Step 1 (preferred): frame index from the gwportal DB instead of a disk scan.

gwportal ingests every py7dt v2 single (ProcessedScienceFrame) with its
path and header metadata, so the frame list comes from the API:

    GET <GWPORTAL_BASE_URL>/api/frames/processed/
        ?version=v2&data_stream=survey&tile_name=T02252,T02385,...
        &filter_name=m400,m425,...&page_size=5000      (X-API-Key header)

Writes the SAME CSV schema as frame_index.py, so select_link.py works
unchanged. Field mapping (gwportal serializer -> index):

    filepath -> path        raw_frame tile -> tile   filter -> band
    unit -> unit            night -> night            mjd -> mjd
    exptime -> exptime      zp (= ZP_AUTO) -> zp      seeing -> seeing
    ul5 (= UL5_AUTO) -> depth, but only if > 0 (UL5_AUTO is 0.0 in v2 headers)
    skyval -> sky, skysig -> skysig, ra/dec_center -> ra_cen/dec_cen
    reduction_status astrom|photom -> has_wcs = 1

A frame whose path is not under one of --mounted (default /lyman; skynet
does not mount /balmer) gets error "not_mounted:<root>", so select_link.py
drops it as unreadable and the count is visible. --check-exists also
stat()s every path.

    export GWPORTAL_BASE_URL=https://proton.snu.ac.kr:8081
    export GWPORTAL_API_KEY=...
    python3 hips7dt/frame_query.py --tiles config/dp0_tiles.txt \
        --out index/dp0_frames.csv
    python3 hips7dt/frame_query.py --demo     # no network
"""

from __future__ import annotations

import argparse
import json
import os
import ssl
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import frame_index as fi  # noqa: E402

MEDIUM_BANDS = [f"m{w}" for w in range(400, 900, 25)]


def http_fetch(base_url: str, api_key: str, params: dict, retries=4) -> dict:
    url = base_url.rstrip("/") + "/api/frames/processed/?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"X-API-Key": api_key,
                                               "Accept": "application/json"})
    ctx = ssl.create_default_context()
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120, context=ctx) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            body = e.read()[:300].decode(errors="replace")
            if e.code < 500 or attempt == retries - 1:
                raise RuntimeError(f"HTTP {e.code} for {url}: {body}") from None
        except urllib.error.URLError:
            if attempt == retries - 1:
                raise
        time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def fetch_all(fetch, params: dict) -> list[dict]:
    """Follow next_page until it is null."""
    out, page = [], 1
    while page:
        data = fetch({**params, "page": page})
        out.extend(data.get("results", []))
        print(f"  page {page}/{data.get('num_pages')}  rows {len(out)}/{data.get('count')}",
              flush=True)
        page = data.get("next_page")
    return out


def to_row(rec: dict, mounted: list[str], check_exists: bool) -> dict:
    row = {f: "" for f in fi.FIELDS}
    path = rec.get("filepath") or ""
    depth = rec.get("ul5")
    row.update({
        "path": path,
        "tile": rec.get("tile") or "",
        "band": rec.get("filter") or "",
        "unit": rec.get("unit") or "",
        "night": rec.get("night") or "",
        "mjd": rec.get("mjd") if rec.get("mjd") is not None else rec.get("raw_frame_mjd"),
        "date_obs": rec.get("obstime") or "",
        "exptime": rec.get("exptime"),
        "zp": rec.get("zp"),
        "seeing": rec.get("seeing"),
        "depth": depth if depth is not None and depth > 0 else None,
        "sky": rec.get("skyval"),
        "skysig": rec.get("skysig"),
        "ra_cen": rec.get("ra_center"),
        "dec_cen": rec.get("dec_center"),
        "has_wcs": int(rec.get("reduction_status") in ("astrom", "photom")),
        "zp_key": "db:zp(ZP_AUTO)", "seeing_key": "db:seeing",
        "depth_key": "db:ul5(UL5_AUTO)", "sky_key": "db:skyval",
    })
    if not path:
        row["error"] = "no_path"
    elif path.endswith("_weight.fits"):
        row["error"] = "weight_map"
    elif mounted and not any(path.startswith(m.rstrip("/") + "/") for m in mounted):
        row["error"] = "not_mounted:" + "/".join(path.split("/")[:2])
    elif check_exists and not os.path.exists(path):
        row["error"] = "missing_on_disk"
    return row


def report(rows):
    import collections
    err = collections.Counter(r["error"] or "ok" for r in rows)
    print(f"frames: {len(rows)}  status: {dict(err)}")
    ok = [r for r in rows if not r["error"] and r["has_wcs"] == 1]
    tb = collections.Counter((r["tile"], r["band"]) for r in ok)
    tiles = sorted({t for t, _ in tb})
    bands = sorted({b for _, b in tb}, key=lambda b: (len(b), b))
    if tiles:
        print("usable frames (readable here, WCS done) per tile x band:")
        print("        " + " ".join(f"{b:>5}" for b in bands))
        for t in tiles:
            print(f"{t:7s} " + " ".join(f"{tb.get((t, b), 0):5d}" for b in bands))
    lost = collections.Counter((r["tile"], r["band"]) for r in rows
                               if r["error"].startswith("not_mounted"))
    if lost:
        print(f"WARNING: {sum(lost.values())} frames are on an unmounted disk "
              f"(e.g. /balmer) in {len(lost)} tile-bands")


def run(a, fetch) -> int:
    tiles = fi.read_list(a.tiles)
    params = {"version": "v2", "data_stream": a.data_stream,
              "filter_name": ",".join(a.bands), "page_size": a.page_size}
    if a.reduction_status:
        params["reduction_status"] = a.reduction_status
    recs = []
    for i in range(0, len(tiles), a.tiles_per_query):
        chunk = tiles[i:i + a.tiles_per_query]
        print(f"query tiles {','.join(chunk)}")
        recs.extend(fetch_all(fetch, {**params, "tile_name": ",".join(chunk)}))
    if a.raw_json:
        Path(a.raw_json).parent.mkdir(parents=True, exist_ok=True)
        Path(a.raw_json).write_text(json.dumps(recs, default=str))
    seen, rows = set(), []
    for rec in recs:
        if rec.get("filepath") in seen:
            continue
        seen.add(rec.get("filepath"))
        rows.append(to_row(rec, a.mounted, a.check_exists))
    fi.write_csv(rows, Path(a.out))
    report(rows)
    print(f"written: {a.out}")
    return 0 if rows else 1


def demo() -> int:
    """Fake API pages shaped like gwportal's _serialize_processed output."""
    def rec(i, tile, band, root="/lyman/data2", status="photom", ul5=0.0):
        return {"id": i, "filename": f"{tile}_{band}_7DT0{i % 3 + 1}_20251213_0741{i:02d}_100s.fits",
                "filepath": f"{root}/processed/2025-12-13/{tile}/{band}/singles/"
                            f"{tile}_{band}_7DT0{i % 3 + 1}_20251213_0741{i:02d}_100s.fits",
                "ra_center": 147.06, "dec_center": -27.1, "mjd": 61022.32 + i * 1e-3,
                "zp": 24.3, "ul5": ul5, "seeing": 1.88, "skyval": 17.0, "skysig": 9.98,
                "tile": tile, "filter": band, "unit": f"7DT0{i % 3 + 1}",
                "night": "2025-12-12", "exptime": 100.0, "obstime": "2025-12-13T07:41:31",
                "reduction_status": status, "processing_version": "v2",
                "data_stream": "survey"}
    recs = [rec(i, t, b) for i, (t, b) in enumerate(
        [(t, b) for t in ("T02524", "T02666") for b in ("m400", "m650")] * 3)]
    recs.append(rec(99, "T02524", "m650", root="/balmer/data1"))
    recs.append(rec(98, "T02666", "m400", status="calib", ul5=17.9))
    recs.append(dict(recs[0]))                            # duplicate row
    pages = [recs[i:i + 5] for i in range(0, len(recs), 5)]
    calls = []

    def fake_fetch(params):
        calls.append(params)
        p = params["page"]
        return {"count": len(recs), "num_pages": len(pages),
                "next_page": p + 1 if p < len(pages) else None,
                "results": pages[p - 1]}

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "idx.csv"
        a = argparse.Namespace(tiles="T02524,T02666", bands=MEDIUM_BANDS,
                               data_stream="survey", reduction_status="",
                               page_size=5000, tiles_per_query=10, raw_json="",
                               out=str(out), mounted=["/lyman"], check_exists=False)
        run(a, fake_fetch)
        import pandas as pd
        df = pd.read_csv(out, dtype=str).fillna("")
        checks = {
            "all pages followed, duplicates dropped": len(calls) == len(pages)
                and len(df) == len(recs) - 1,
            "query asks for v2 survey, tiles and bands": calls[0]["version"] == "v2"
                and calls[0]["tile_name"] == "T02524,T02666"
                and "m875" in calls[0]["filter_name"],
            "balmer frame flagged not_mounted": (df["error"] == "not_mounted:/balmer").sum() == 1,
            "calib-only frame has no WCS": df.loc[df["path"].str.contains("_074198_"), "has_wcs"].tolist() == ["0"],
            "UL5_AUTO 0.0 -> depth empty": (df.loc[df["has_wcs"] == "1", "depth"] == "").all(),
        }
        # the index feeds select_link.py unchanged
        import select_link as sl
        sa = argparse.Namespace(index=str(out), out=str(Path(tmp) / "tree"), bands=MEDIUM_BANDS,
                                tiles=None, exclude_dates="", exclude_images=[],
                                require_zp=True, max_seeing=0, per_tile_band=10,
                                zp_warn=0.1, scaled_copy=None, dry_run=True,
                                zp_clip=0.3, hhh_scale=None)
        sl.run(sa)
        man = pd.read_csv(Path(tmp) / "tree/manifest.csv", dtype=str)
        checks["select_link reads it (12 selected)"] = (man["selected"] == "1").sum() == 12
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles", help="file or comma list of tiles")
    ap.add_argument("--bands", default=",".join(MEDIUM_BANDS))
    ap.add_argument("--out", default="index/frames.csv")
    ap.add_argument("--data-stream", default="survey")
    ap.add_argument("--reduction-status", default="",
                    help="e.g. photom or astrom,photom; default all (WCS judged from it)")
    ap.add_argument("--mounted", default="/lyman",
                    help="comma list of path roots readable here; '' = no check")
    ap.add_argument("--check-exists", action="store_true")
    ap.add_argument("--page-size", type=int, default=5000)
    ap.add_argument("--tiles-per-query", type=int, default=10)
    ap.add_argument("--raw-json", default="", help="also save the raw API records")
    ap.add_argument("--base-url", default=os.getenv("GWPORTAL_BASE_URL", ""))
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        return demo()
    if not a.tiles:
        ap.error("--tiles is required")
    key = os.getenv("GWPORTAL_API_KEY")
    if not (a.base_url and key):
        ap.error("set GWPORTAL_BASE_URL (or --base-url) and GWPORTAL_API_KEY")
    a.bands = [b for b in a.bands.replace(",", " ").split() if b]
    a.mounted = [m for m in a.mounted.split(",") if m]
    return run(a, lambda p: http_fetch(a.base_url, key, p))


if __name__ == "__main__":
    sys.exit(main())
