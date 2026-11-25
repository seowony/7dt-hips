# CLAUDE.md

Working notes for Claude Code sessions on this repository. Code and repo
documents are English; discussion with the maintainer is Korean, replies
concise.

## What this is

Colour HiPS of 7DT medium-band data. Per-band FITS HiPS are built with
Hipsgen from **single processed science frames of the v2 reprocessing**
(not the coadds), then combined into a colour PNG HiPS with the recipe of
Thomas Boch (CDS). Thomas's test result is at
https://alasky.cds.unistra.fr/7DT/test/ ; his script is kept verbatim in
`docs/thomas/make_color_hips_thomas_v1.py`.

Decisions (maintainer, 2026-11-25):
- Thomas's colour parameters are the **default**; change them only after
  the default result has been seen, and record why.
- Scope = the RIS DP0/DP1/DP2 footprints. **DP0 (7 IMS tiles) first**, as
  the example run.
- Input frames are used through **symlinks** (no copies).
- The finished HiPS will be **linked from gwportal** (maintainer,
  2026-11-25), so the output must end up where gwportal's web server can
  serve it.
- DP0 uses the **best N frames per tile x band** (by seeing), not all of
  them (~600 per tile-band exist): output size does not depend on N, but
  Hipsgen time does, and DP0 should not be ~3 mag deeper than DP1/DP2.
  N = 10 vs 30 is compared on m650 first.
- **All data products live under `/data8/skyteam/7DT`** on skynet
  (maintainer, 2026-11-25): `index/`, `hips_in/`, `hips/`, `logs/`,
  `hipsgen_cache/`. `hips7dt/paths.py`: relative paths given to any script
  (and IN/OUT of `run_hipsgen.sh`) resolve under `$HIPS7DT_WORK` (default
  `/data8/skyteam/7DT`), absolute paths are used as is. The git checkout
  (`~/projects/7dt-hips`) holds code only.
- Runs on **skynet**: `/lyman/data1`, `/lyman/data2` are mounted there;
  `/balmer` is NOT. Frames that live only on balmer are invisible from
  skynet - check per tile x band counts in the index.

## Where Claude runs vs. where scripts run

A Claude Code cloud session cannot reach skynet, proton or /lyman. Claude
writes code and validates it with each script's `--demo` (synthetic data,
no server); the maintainer runs it on skynet and pastes back logs /
reports. **Every new script must ship with a `--demo` self-check.**
`make_color_hips.py --demo` checks pixel identity with Thomas's original;
keep that check passing unless the recipe is changed on purpose.

## Pipeline

1. `hips7dt/frame_query.py` - frame index from the **gwportal DB** via
   `GET /api/frames/processed/?version=v2&data_stream=survey&tile_name=...`
   (maintainer, 2026-11-25: frame paths are taken from the DB that gwportal
   serves, not from a disk scan). Same CSV schema as `frame_index.py` (the
   header-scan fallback). Frames not under `/lyman` are flagged
   `not_mounted`. Mapping source: gwportal `api/views.py::_serialize_processed`
   and `survey/models.py::ProcessedScienceFrame` (zp = ZP_AUTO, ul5 = UL5_AUTO).
2. `hips7dt/select_link.py` - per tile x band best-N selection, exclusions,
   symlink tree, `manifest.csv`, `zp_report.csv`; `--scaled-copy ZPREF`
   when ZPs differ (Hipsgen cannot scale individual files).
3. `hips7dt/run_hipsgen.sh` - one FITS HiPS per band, order 10 (0.40"/pix).
4. `hips7dt/make_color_hips.py` - colour PNG HiPS (resume by default).

## Open items

- [x] v2 frame path pattern and header keywords - from gwportal
      (`ingest_processed_frames_v2.py`, `FITS_HEADERS/v2_processedsingle.head`).
- [ ] Is the v2 single background-subtracted? (SKYVAL ~17 ADU suggests yes.)
- [ ] Does Hipsgen read the TPV WCS of v2 singles correctly?
- [ ] gwportal's API does not return `is_production` / `sanity_flag`; add them
      to its serializer if QA flags are needed for selection.
- [x] Hipsgen parameter names checked against `-h` (v12.903): note
      `skyVal` (capital V), `tileWidth`, `expTime`, `mode`.
- [ ] Hipsgen follows file symlinks? (first m650 run)
- [ ] ZP spread. First DP0 selection (N=10, 1,400 frames): per-band ZP std
      0.03-0.13 mag, but outliers up to 7 mag (m625 22.8-29.9, m675 21.0,
      m850 28.1, m450/m600 ~25.05 vs ~26.5) -> `--zp-clip 0.3` (default)
      drops them before ranking. Remaining offsets: try `--hhh-scale ZPREF`
      (`<name>.hhh` header with HIPSSCAL beside each symlink) + Hipsgen
      `expTime=HIPSSCAL` (SCALE=1) - no pixel copies. Needs confirming that
      Hipsgen reads .hhh headers and divides by expTime (`-man expTime`).
- [ ] DP0 audit problem images (2025-06/07 m425/m625/m850, 2025-09-04 all
      bands, 2025-10 red bands) were found on coadds; check whether v2
      frames share them. Lists on proton: ris-db-pipeline
      `outputs/release_integrity/{bands_dp0_culprits,epoch_bad_images}.csv`.
- [ ] After DP0: DP1/DP2 footprints (~1,500 deg^2, ~13 TB of order-10 FITS
      HiPS for 20 bands) - process the sky in chunks and delete the FITS
      after the PNG, or keep FITS at order 9.

## Workflow

Claude works on a session branch, pushes, and opens a pull request; the
maintainer merges it on GitHub. Always end the reply with the PR's merge
page link (https://github.com/seowony/7dt-hips/pull/N). Delete merged
branches. Write decisions made in chat into this file before a session
ends - it is the only context a new session inherits.
