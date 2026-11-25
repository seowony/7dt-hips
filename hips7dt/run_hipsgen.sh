#!/usr/bin/env bash
# Step 3: one FITS HiPS per band with Hipsgen, from the symlink tree of
# select_link.py. Relative IN/OUT are taken under $HIPS7DT_WORK
# (default /data8/skyteam/7DT). Output layout is what make_color_hips.py expects:
#   $OUT/m400/{properties,Moc.fits,Norder*/...}, $OUT/m425/..., ...
#
# Bands run one after another; Hipsgen itself uses $THREADS threads.
# A band is skipped when $OUT/<band>/properties AND Norder3/Allsky.fits exist
# (resume). Wall time, Hipsgen exit code and completeness per band go to
# $OUT/timing.tsv; a non-zero exit code does not stop the other bands.
#
# Parameter names checked against `Hipsgen.jar -h` (Aladin v12.903):
#   in= out= id= order= tileWidth= bitpix= maxThread= skyVal= expTime= mode=
#
#   IN=hips_in/dp0_n10 OUT=hips/dp0_n10 BANDS=m650 bash hips7dt/run_hipsgen.sh
#   ZP / sky: do them in select_link.py (--scaled-copy 25 --sky-sub). Tested
#   2026-11-25 on m650: Hipsgen expTime= with a .hhh header changed nothing,
#   skyVal=auto multiplied the pixels by ~2.25 instead of removing the sky.
#   SKY (default SKYVAL) -> skyVal=SKYVAL: Hipsgen subtracts each frame's own
#   header SKYVAL while coadding (files untouched). Tested 2026-11-25, m650
#   T02524 region: sky 58 -> 1.4 ADU, star excess unchanged. SKY= (empty)
#   switches it off.
#   EXTRA="mode=..."  any further Hipsgen parameters
set -euo pipefail

HIPSGEN=${HIPSGEN:-$HOME/bin/Hipsgen.jar}
WORK=${HIPS7DT_WORK:-/data8/skyteam/7DT}   # relative IN/OUT are under WORK
IN=${IN:?set IN=<symlink tree from select_link.py>}
OUT=${OUT:?set OUT=<HiPS output root>}
case "$IN" in /*) ;; *) IN="$WORK/$IN" ;; esac
case "$OUT" in /*) ;; *) OUT="$WORK/$OUT" ;; esac
LOGS=${LOGS:-$WORK/logs}
ORDER=${ORDER:-10}                 # 512-px tiles at order 10 -> 0.40"/pix (7DT 0.505")
THREADS=${THREADS:-8}               # maintainer: use at most 8 cores on skynet
MEM=${MEM:-16g}
ID_PREFIX=${ID_PREFIX:-SNU/P/7DT/RIS-DP0}
BANDS=${BANDS:-$(ls "$IN" | grep -E '^m[0-9]{3}$')}
SKY=${SKY-SKYVAL}
EXTRA=${EXTRA:-}

opts=(order="$ORDER" tileWidth=512 bitpix=-32 maxThread="$THREADS" cache="$WORK/hipsgen_cache")
[ -n "$SKY" ] && opts+=(skyVal="$SKY")
# shellcheck disable=SC2206
[ -n "$EXTRA" ] && opts+=($EXTRA)

mkdir -p "$OUT" "$LOGS"
echo "IN=$IN  OUT=$OUT  LOGS=$LOGS"
for band in $BANDS; do
  [ -d "$IN/$band" ] || { echo "skip $band: no $IN/$band"; continue; }
  if [ -f "$OUT/$band/properties" ] && [ -f "$OUT/$band/Norder3/Allsky.fits" ]; then
    echo "done already: $band"; continue
  fi
  n=$(find -L "$IN/$band" -name '*.fits' | wc -l)
  echo "=== $band: $n frames  ${opts[*]}  $(date)"
  t0=$(date +%s)
  # Hipsgen can exit non-zero after a complete run (e.g. warnings); with
  # `set -e -o pipefail` that used to stop the whole loop after the first
  # band. Record the exit code instead and judge completion by the files.
  set +e
  java -Xmx"$MEM" -jar "$HIPSGEN" \
      in="$IN/$band" out="$OUT/$band" id="$ID_PREFIX/$band" "${opts[@]}" \
      2>&1 | tee "$LOGS/hipsgen_$(basename "$OUT")_$band.log"
  rc=${PIPESTATUS[0]}
  set -e
  done_ok=no
  [ -f "$OUT/$band/properties" ] && [ -f "$OUT/$band/Norder3/Allsky.fits" ] && done_ok=yes
  printf '%s\t%s\t%d\trc=%s\tcomplete=%s\t%s\n' "$band" "$n" "$(( $(date +%s) - t0 ))" \
      "$rc" "$done_ok" "${opts[*]}" >> "$OUT/timing.tsv"
  echo "=== $band finished: rc=$rc complete=$done_ok"
done
