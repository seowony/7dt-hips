#!/usr/bin/env bash
# Step 3: one FITS HiPS per band with Hipsgen, from the symlink tree of
# select_link.py. Output layout is what make_color_hips.py expects:
#   $OUT/m400/{properties,Moc.fits,Norder*/...}, $OUT/m425/..., ...
#
# Bands run one after another; Hipsgen itself uses $THREADS threads.
# A band whose $OUT/<band>/properties exists is skipped (resume; needs properties AND Norder3/Allsky.fits).
#
# Parameter names follow Hipsgen's usual syntax but were NOT checked against
# the installed version - run `java -jar $HIPSGEN -h` once and adjust
# EXTRA / the line below if a name differs.
#
#   IN=hips_in/dp0 OUT=hips/dp0 bash hips7dt/run_hipsgen.sh            # all bands
#   BANDS="m650" IN=hips_in/dp0 OUT=hips/dp0 bash hips7dt/run_hipsgen.sh  # one band
set -euo pipefail

HIPSGEN=${HIPSGEN:-$HOME/bin/Hipsgen.jar}
IN=${IN:?set IN=<symlink tree from select_link.py>}
OUT=${OUT:?set OUT=<HiPS output root>}
ORDER=${ORDER:-10}                 # 512-px tiles at order 10 -> 0.40"/pix (7DT 0.505")
THREADS=${THREADS:-4}
MEM=${MEM:-16g}
ID_PREFIX=${ID_PREFIX:-SNU/P/7DT/RIS-DP0}
BANDS=${BANDS:-$(ls "$IN")}
EXTRA=${EXTRA:-"skyval=auto"}     # per-frame sky subtraction; set EXTRA="" to compare

mkdir -p "$OUT" logs
for band in $BANDS; do
  [ -d "$IN/$band" ] || { echo "skip $band: no $IN/$band"; continue; }
  if [ -f "$OUT/$band/properties" ] && [ -f "$OUT/$band/Norder3/Allsky.fits" ]; then
    echo "done already: $band"; continue
  fi
  n=$(find -L "$IN/$band" -name '*.fits*' | wc -l)
  echo "=== $band: $n frames  $(date)"
  java -Xmx"$MEM" -jar "$HIPSGEN" \
      in="$IN/$band" out="$OUT/$band" id="$ID_PREFIX/$band" \
      order="$ORDER" bitpix=-32 maxThread="$THREADS" $EXTRA \
      2>&1 | tee "logs/hipsgen_$(basename "$OUT")_$band.log"
done
