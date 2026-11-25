# 7dt-hips — 7DT colour HiPS (RIS DP0 first)

Builds per-band FITS HiPS from 7DT single processed science frames (v2
reprocessing) and a colour PNG HiPS with the colour recipe of Thomas Boch
(CDS), whose parameter values are kept as the defaults. The first test
covers the **seven RIS DP0 IMS tiles** and runs on **skynet** (`/lyman` is
mounted there, `/balmer` is not).

## Steps

| step | script | input → output |
|---|---|---|
| 1 | `hips7dt/frame_query.py` | gwportal DB (API) → `index/dp0_frames.csv` |
| 1′ | `hips7dt/frame_index.py` | fallback: disk scan of frame headers → same CSV |
| 2 | `hips7dt/select_link.py` | index → `hips_in/dp0/<band>/` symlink tree + `manifest.csv`, `zp_report.csv` |
| 3 | `hips7dt/run_hipsgen.sh` | symlink tree → `hips/dp0/<band>/` FITS HiPS (Hipsgen) |
| 4 | `hips7dt/make_color_hips.py` | `hips/dp0/m*/` → `hips/dp0/color/` PNG HiPS |

Every Python step has `--demo`, a self-check that needs no data.
`make_color_hips.py --demo` also checks that it writes **pixel-identical
tiles** to Thomas's original script (`docs/thomas/`) from the same input.

## Where the data go

Every data product is written under `$HIPS7DT_WORK` (default
`/data8/skyteam/7DT`): `index/`, `hips_in/` (symlink trees), `hips/` (HiPS),
`logs/`, `hipsgen_cache/`. Relative paths given to the scripts are taken
under it; absolute paths are used as is. Run the scripts from the git
checkout, which holds code only.

## Setup on skynet (once)

```bash
pip install astropy numpy pandas mocpy pillow
mkdir -p ~/bin && wget -O ~/bin/Hipsgen.jar https://aladin.cds.unistra.fr/java/Hipsgen.jar
java -jar ~/bin/Hipsgen.jar -h | less     # check the parameter names (see Caveats)
python3 hips7dt/frame_index.py --demo && python3 hips7dt/select_link.py --demo \
  && python3 hips7dt/make_color_hips.py --demo
```

## DP0 run

```bash
# 1) index from the gwportal DB (v2 survey singles of the DP0 tiles)
export GWPORTAL_BASE_URL=https://proton.snu.ac.kr:8081
export GWPORTAL_API_KEY=...            # same key as toolkit/gwportal_client.py
python3 hips7dt/frame_query.py --tiles config/dp0_tiles.txt \
    --out index/dp0_frames.csv --check-exists
#    prints frames per tile x band and how many sit on /balmer (not on skynet)

# 2) selection: best N frames by seeing per tile x band, ZP outliers dropped.
#    Hipsgen cannot rescale single files, so write float32 copies brought to a
#    common ZP (25) with each frame's SKYVAL subtracted (~245 MB per frame).
python3 hips7dt/select_link.py --index index/dp0_frames.csv --out hips_in/dp0_n10_zp \
    --tiles config/dp0_tiles.txt --exclude-dates 2025-09-04 --per-tile-band 10 \
    --require-zp --scaled-copy 25 --sky-sub --workers 8

# 3) Hipsgen. Test on a small region first (EXTRA="region=<moc.fits>"),
#    then all bands; THREADS up to what the machine allows
THREADS=16 IN=hips_in/dp0_n10_zp OUT=hips/dp0_n10_zp bash hips7dt/run_hipsgen.sh

# 4) colour composite
python3 hips7dt/make_color_hips.py --in-root hips/dp0_n10 --out hips/dp0_n10/color --workers 4
```

To inspect: `python3 hips7dt/make_viewer.py` writes `hips/index.html`, two
linked Aladin Lite panels with a drop-down of every HiPS under `hips/`.
Serve it on skynet (`cd /data8/skyteam/7DT/hips && python3 -m http.server
8090 --bind 127.0.0.1`), tunnel from the laptop (`ssh -N -L
8090:127.0.0.1:8090 skyteam@skynet`) and open http://localhost:8090/.
Aladin Desktop can open the same URLs (File > Open > URL).

## v2 single frames (from gwportal)

- Path: `/lyman/data2/processed/{date}/{tile}/{filter}/singles/
  {tile}_{filter}_{unit}_{date}_{time}_{exp}s.fits`, with a
  `*_weight.fits` weight map beside it (skipped). The same layout exists on
  `/balmer/data1/processed/...`, which skynet cannot read.
- Header sample: gwportal `FITS_HEADERS/v2_processedsingle.head`. 9576 x 6388
  px, TPV WCS, `ZP_AUTO` (DB `zp`), `SEEING`, `SKYVAL`/`SKYSIG`, `UL5_0..5`;
  `UL5_AUTO` is 0.0 and is treated as missing.
- `SKYVAL` ~17 ADU in a 100-s frame suggests the background is already
  subtracted - to be confirmed before choosing Hipsgen's `skyval`.

## Caveats

- **Hipsgen parameters** (`in= out= id= order= bitpix= maxThread=
  skyval=auto`) have not been run yet. If they differ from the `-h`
  output, fix `run_hipsgen.sh`. The colour script needs `hips_pixel_cut`
  in each band's `properties`, plus `Norder3/Allsky.fits` and `Moc.fits`.
- **Symlinks**: check with one band (m650) that Hipsgen follows *file*
  symlinks (number of files read in the log).
- **ZP spread**: Hipsgen cannot scale individual files. If a band's
  `zp_range` in `zp_report.csv` exceeds ~0.1 mag (unit-to-unit
  sensitivity), the mosaic will show seams, so use `--scaled-copy 25
  --sky-sub` (float32 copies; the disk cost is printed first). Hipsgen's own
  `expTime=` (with a .hhh header) and `skyVal=auto` did not do this (tested).
  Also check in the header whether the ZP includes the exposure time.
- **DP0 has many dates** (140–170 per tile). Do not use all of them;
  cut with `--per-tile-band`. 10 frames × 7 tiles × 20 bands ≈ 1,400 frames.
- **Problem images from the DP0 audit** (2025-06/07 m425, m625, m850;
  2025-09-04 all bands; 2025-10 red bands) were found on the coadds. It
  is unknown whether v2 fixes them, so only 2025-09-04 is excluded at
  first; check those bands in the result. The audit lists on proton
  (`outputs/release_integrity/{bands_dp0_culprits,epoch_bad_images}.csv`)
  can be copied over and passed as `--exclude-images a.csv,b.csv`.
- **Disk**: 7 tiles ≈ 8–9 deg²; an order-10 FITS HiPS ≈ 0.43 GB/deg² per
  band → ~4 GB per band, ~80 GB for 20 bands.

## Improvements (after the default result is in)

1. Gradients: compare with and without `skyval=auto`; check whether the
   v2 frames are already background-subtracted.
2. Common ZP: does `--scaled-copy` reduce the seams?
3. Colour fringes where band coverage ends; brightness jumps between
   orders (Thomas's per-order blackpoint).
4. Compare with a HiPS built from the coadds of the same field (cosmic
   rays, satellite trails).
