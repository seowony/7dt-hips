#!/usr/bin/env python3
"""Write an Aladin Lite comparison page for every HiPS under a root.

Finds each directory holding a HiPS `properties` file below --root (default
$HIPS7DT_WORK/hips) and writes <root>/index.html with two linked Aladin
Lite panels (same position / zoom), each with a drop-down of all HiPS, so
two variants (e.g. test_raw/m650 vs test_scaled/m650) can be compared side
by side. Serve the root over HTTP and open the page:

    # skynet
    python3 hips7dt/make_viewer.py
    cd /data8/skyteam/7DT/hips && python3 -m http.server 8090 --bind 127.0.0.1
    # laptop
    ssh -N -L 8090:127.0.0.1:8090 skyteam@skynet
    # browser: http://localhost:8090/

Aladin Lite is loaded from the CDS site, so the laptop needs internet.
The tiles shown are Hipsgen's PNG previews (its own pixel cut).

    python3 hips7dt/make_viewer.py --demo
"""

from __future__ import annotations

import argparse
import html
import json
import sys
import tempfile
from pathlib import Path

PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>7DT HiPS viewer</title>
<script src="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"></script>
<style>
 body {{ margin: 0; font: 14px sans-serif; background: #111; color: #ddd; }}
 header {{ padding: 6px 10px; }}
 .row {{ display: flex; gap: 4px; height: calc(100vh - 80px); }}
 .col {{ flex: 1; display: flex; flex-direction: column; }}
 .col select {{ margin: 2px 0; }}
 .al {{ flex: 1; }}
</style></head><body>
<header>7DT HiPS under <code>{root}</code> - {n} HiPS. Left and right panels move together.
 Target: <input id="tgt" value="{target}" size="22"> <button id="go">go</button></header>
<div class="row">
 <div class="col"><select id="s0"></select><div id="a0" class="al"></div></div>
 <div class="col"><select id="s1"></select><div id="a1" class="al"></div></div>
</div>
<script>
const HIPS = {hips_json};
A.init.then(() => {{
  const opts = {{target: "{target}", fov: {fov}, showReticle: false,
                 showLayersControl: true, showFullscreenControl: false}};
  const al = [A.aladin('#a0', opts), A.aladin('#a1', opts)];
  ['s0', 's1'].forEach((id, i) => {{
    const sel = document.getElementById(id);
    HIPS.forEach((h, k) => sel.add(new Option(h, h, false, k === Math.min(i, HIPS.length - 1))));
    const show = () => al[i].setImageSurvey(new URL(sel.value + '/', location.href).href);
    sel.onchange = show; show();
  }});
  let busy = false;
  const sync = (src, dst) => () => {{
    if (busy) return; busy = true;
    const [ra, dec] = al[src].getRaDec(); al[dst].gotoRaDec(ra, dec);
    al[dst].setFoV(al[src].getFov()[0]); busy = false;
  }};
  al[0].on('positionChanged', sync(0, 1)); al[0].on('zoomChanged', sync(0, 1));
  al[1].on('positionChanged', sync(1, 0)); al[1].on('zoomChanged', sync(1, 0));
  document.getElementById('go').onclick = () => al[0].gotoObject(document.getElementById('tgt').value);
}});
</script></body></html>
"""


def find_hips(root: Path) -> list[str]:
    out = []
    for prop in sorted(root.rglob("properties")):
        d = prop.parent
        if d.name == "HpxFinder":   # Hipsgen's index, not an image HiPS
            continue
        if (d / "Norder3").exists() or any(d.glob("Norder*")):
            out.append(d.relative_to(root).as_posix())
    return out


def write_page(root: Path, target: str, fov: float) -> Path:
    hips = find_hips(root)
    page = root / "index.html"
    page.write_text(PAGE.format(root=html.escape(str(root)), n=len(hips),
                                target=html.escape(target), fov=fov,
                                hips_json=json.dumps(hips)))
    return page


def demo() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for name in ("test_raw/m650", "test_scaled/m650", "dp0_n10/color"):
            (root / name / "Norder3").mkdir(parents=True)
            (root / name / "properties").write_text("hips_order = 10\n")
        (root / "not_hips").mkdir()
        (root / "test_raw/m650/HpxFinder/Norder3").mkdir(parents=True)
        (root / "test_raw/m650/HpxFinder/properties").write_text("x = 1\n")
        page = write_page(root, "79.151 -60.287", 0.5)
        text = page.read_text()
        checks = {
            "three HiPS found, others ignored": find_hips(root) == [
                "dp0_n10/color", "test_raw/m650", "test_scaled/m650"],
            "page lists them and loads Aladin Lite": all(
                h in text for h in ("test_raw/m650", "aladin.js", "79.151 -60.287")),
        }
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if all(checks.values()) else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="hips")
    ap.add_argument("--target", default="79.151 -60.287", help="initial 'ra dec' or name")
    ap.add_argument("--fov", type=float, default=0.5, help="initial field of view [deg]")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo:
        return demo()
    sys.path.insert(0, str(Path(__file__).parent))
    from paths import resolve
    root = Path(resolve(a.root))
    page = write_page(root, a.target, a.fov)
    print(f"{len(find_hips(root))} HiPS -> {page}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
