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

1. `hips7dt/frame_index.py` - header-only index of v2 frames (candidate
   keyword lists for ZP / seeing / depth / sky; `--dump-header`).
2. `hips7dt/select_link.py` - per tile x band best-N selection, exclusions,
   symlink tree, `manifest.csv`, `zp_report.csv`; `--scaled-copy ZPREF`
   when ZPs differ (Hipsgen cannot scale individual files).
3. `hips7dt/run_hipsgen.sh` - one FITS HiPS per band, order 10 (0.40"/pix).
4. `hips7dt/make_color_hips.py` - colour PNG HiPS (resume by default).

## Open items

- [ ] v2 frame path pattern and header keywords (ZP, seeing, depth, sky,
      background-subtracted or not) - from `--dump-header` on skynet.
- [ ] Hipsgen parameter names not yet checked against `Hipsgen.jar -h`
      (`in= out= id= order= bitpix= maxThread= skyval=auto`), and whether it
      follows file symlinks.
- [ ] ZP spread per band (`zp_report.csv`) -> symlinks or scaled copies.
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
