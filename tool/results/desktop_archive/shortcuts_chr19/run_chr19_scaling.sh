#!/bin/bash
# MEMs of two chromosome 19 haplotypes (HG00423 and HG00438, haplotype 1)
# against growing collections of chromosome 19s from HPRC release 2: CHM13
# as the reference plus both haplotypes of the first 16, 32 and 64
# individuals by sample ID (each collection a prefix of the next, the
# query individuals left out of all of them). For each size: the RLBWTs of
# the collection and of the reversed collection (Big-BWT), shortcuts for
# ell = 16, 25, 50, tau = 3, 6, theta = 0, and BML (L = 100, 1000) and
# forward-backward (L = 1000), each search timed 5 times with and 5 times
# without shortcuts, alternating, keeping the medians.
#
# Put this script and reportable-shortcuts-code.tar.gz in the same folder
# (e.g. ~/research/rlbwt-tricks/shortcuts), close other programs, including the Claude app, and
# run, with DATADIR a folder on the external drive:
#   nohup bash ~/research/rlbwt-tricks/shortcuts/run_chr19_scaling.sh DATADIR > ~/research/rlbwt-tricks/shortcuts/run_chr19_scaling.log 2>&1 &
# If it is stopped, the same command resumes it: finished downloads,
# indexes and query results are kept.
#
# Environment variables: SIZES=16,32,64 (individuals per collection),
# QUERY_SAMPLES=HG00423,HG00438, REPEAT=5, THREADS=4 (for compiling),
# WORK=~/shortcuts_build (where the code and Big-BWT are built; not on the
# external drive, which may not allow running programs), BIGBWT=<path to
# Big-BWT's bigbwt script> (otherwise ~/Big-BWT/bigbwt, or it is cloned
# from GitHub and built).
#
# Downloads about 4 GB. Disk on DATADIR: about 8 GB for the downloaded
# collection, plus about 3 times the size of the collection being indexed
# while it is indexed (about 25 GB for 64 individuals), plus a few GB of
# indexes per size.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
if [ $# -lt 1 ]; then sed -n '2,30p' "$0"; exit 1; fi
mkdir -p "$1"
DATA=$(cd "$1" && pwd)
SIZES=${SIZES:-16,32,64}
QUERY_SAMPLES=${QUERY_SAMPLES:-HG00423,HG00438}
REPEAT=${REPEAT:-5}
THREADS=${THREADS:-4}
WORK=${WORK:-$HOME/shortcuts_build}
PREPARE_ARGS=${PREPARE_ARGS:-}  # for testing
stamp() { echo "[$(date '+%F %T')] $*"; }
free_gb() { df -BG --output=avail "$1" | tail -1 | tr -dc 0-9; }
fail() { echo "error: $*" >&2; exit 1; }

for tool in python3 cmake g++ make taskset; do
  command -v "$tool" > /dev/null || fail "$tool is not installed (sudo apt install cmake g++ make python3 util-linux)"
done
[ -r "$HERE/reportable-shortcuts-code.tar.gz" ] || fail "put reportable-shortcuts-code.tar.gz next to $0"
MAX=0
for s in ${SIZES//,/ }; do [ "$s" -gt "$MAX" ] && MAX=$s; done
stamp "sizes $SIZES, queries $QUERY_SAMPLES, data in $DATA ($(free_gb "$DATA") GB free)"
free -g | head -2

stamp "building and testing the code"
mkdir -p "$WORK"
rm -rf "$WORK/reportable-shortcuts"
tar xzf "$HERE/reportable-shortcuts-code.tar.gz" -C "$WORK"
cmake -S "$WORK/reportable-shortcuts" -B "$WORK/reportable-shortcuts/build" -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build "$WORK/reportable-shortcuts/build" -j"$THREADS" > /dev/null
B="$WORK/reportable-shortcuts/build"
S="$WORK/reportable-shortcuts/scripts"
"$B/test_bitvector"
"$B/test_shortcuts" 200 | tail -1
"$B/test_mems" 200 | tail -1

if [ -z "${BIGBWT:-}" ]; then
  if [ -x "$HOME/Big-BWT/bigbwt" ]; then
    BIGBWT="$HOME/Big-BWT/bigbwt"
  else
    if [ ! -x "$WORK/Big-BWT/bigbwt" ]; then
      stamp "building Big-BWT"
      command -v git > /dev/null || fail "git is not installed (sudo apt install git), or set BIGBWT"
      rm -rf "$WORK/Big-BWT"
      git clone -q https://github.com/alshai/Big-BWT.git "$WORK/Big-BWT"
      make -C "$WORK/Big-BWT" > /dev/null 2>&1 || fail "building Big-BWT failed (it needs zlib: sudo apt install zlib1g-dev)"
    fi
    BIGBWT="$WORK/Big-BWT/bigbwt"
  fi
fi
stamp "Big-BWT is $BIGBWT"
TIME=()
if [ -x /usr/bin/time ]; then TIME=(/usr/bin/time -v); fi

D="$DATA/hprc_chr19"
if [ ! -s "$D/collection.txt" ]; then
  stamp "downloading chromosome 19s for $MAX individuals"
  # shellcheck disable=SC2086
  python3 "$S/prepare_hprc_chr19.py" "$D" --individuals "$MAX" --query-samples "$QUERY_SAMPLES" $PREPARE_ARGS
fi
cut -f1-7 "$D/manifest.tsv"
N=$(cat "$D/ref_length.txt")

# The CPU with the highest maximum frequency.
PIN=$(lscpu -p=CPU,MAXMHZ 2>/dev/null | grep -v '^#' | sort -t, -k2 -gr | head -1 | cut -d, -f1 || true)
PIN=${PIN:-0}
stamp "queries pinned to CPU $PIN"

# A directory of shortcuts is complete if it has all three files and three
# lines in build.csv after the header; a CSV of MEMs if it has a row per
# pattern and minimum length.
complete() { [ -s "$1/rlbwt.bin" ] && [ -s "$1/build.csv" ] && [ "$(wc -l < "$1/build.csv")" -ge 4 ] &&
             [ -s "$1/ell16_theta0.shortcuts" ] && [ -s "$1/ell25_theta0.shortcuts" ] && [ -s "$1/ell50_theta0.shortcuts" ]; }
rows() { [ -s "$1" ] && echo $(( $(wc -l < "$1") - 1 )) || echo 0; }
QUERIES=$(grep -c . "$D/queries.txt")

for size in ${SIZES//,/ }; do
  X="$D/n$size"
  mkdir -p "$X/mems"
  stamp "=== $size individuals ($(free_gb "$DATA") GB free)"
  if ! complete "$X/fwd_t3" || ! complete "$X/fwd_t6" || ! complete "$X/rev_t3" || ! complete "$X/rev_t6"; then
    # The collection: the reference and 2 * size haplotypes.
    if [ "$size" -eq "$MAX" ]; then
      C="$D/collection.txt"
    else
      C="$X/collection.txt"
      head -n $((1 + 2 * size)) "$D/collection.txt" > "$C"
    fi
    gb=$(( $(stat -c %s "$C") / 1000000000 + 1 ))
    [ "$(free_gb "$DATA")" -ge $((3 * gb + 5)) ] || fail "need about $((3 * gb + 5)) GB free on $DATA to index $size individuals"
    for side in fwd rev; do
      if complete "$X/${side}_t3" && complete "$X/${side}_t6"; then continue; fi
      if [ "$side" = fwd ]; then IN="$C"; FLAG=(); else IN="$X/collection_rev.txt"; FLAG=(--reversed); fi
      if ! complete "$X/${side}_t3"; then
        rm -rf "$X/${side}_t3"
        if [ "$side" = rev ]; then python3 "$S/reverse_collection.py" "$C" "$IN"; fi
        stamp "Big-BWT, $side"
        rm -f "$IN".*
        ${TIME[@]+"${TIME[@]}"} "$BIGBWT" "$IN" || true
        [ -s "$IN.bwt" ] || { tail -5 "$IN.log" >&2; fail "Big-BWT failed on $IN"; }
        stamp "RLBWT and shortcuts, $side, tau = 3"
        "$B/build_index" --text "$C" --bwt "$IN.bwt" --reference-length "$N" --out "$X/${side}_t3" \
          --ell 16,25,50 --theta 0 --max-tier 3 ${FLAG[@]+"${FLAG[@]}"}
        rm -f "$IN".*
        if [ "$side" = rev ]; then rm -f "$IN"; fi
      fi
      stamp "shortcuts, $side, tau = 6"
      rm -rf "$X/${side}_t6"
      "$B/build_index" --text "$C" --rlbwt "$X/${side}_t3/rlbwt.bin" --reference-length "$N" --out "$X/${side}_t6" \
        --ell 16,25,50 --theta 0 --max-tier 6 ${FLAG[@]+"${FLAG[@]}"}
    done
    [ "$C" = "$D/collection.txt" ] || rm -f "$C"
  fi
  cat "$X"/fwd_t3/build.csv "$X"/fwd_t6/build.csv "$X"/rev_t3/build.csv "$X"/rev_t6/build.csv

  for t in 6 3; do
    for ell in 25 16 50; do
      c="t${t}_ell${ell}"
      if [ "$(rows "$X/mems/${c}_bml.csv")" -lt $((2 * QUERIES)) ]; then
        stamp "$size individuals, $c: BML"
        taskset -c "$PIN" "$B/query_mems" --index "$X/fwd_t$t" --reversed-index "$X/rev_t$t" --shortcuts "ell${ell}_theta0" \
          --patterns "$D/queries.txt" --algorithm bml --min-length 100,1000 --repeat "$REPEAT" --csv "$X/mems/${c}_bml.csv"
      fi
      if [ "$(rows "$X/mems/${c}_fb.csv")" -lt "$QUERIES" ]; then
        stamp "$size individuals, $c: forward-backward"
        taskset -c "$PIN" "$B/query_mems" --index "$X/fwd_t$t" --reversed-index "$X/rev_t$t" --shortcuts "ell${ell}_theta0" \
          --patterns "$D/queries.txt" --algorithm fb --min-length 1000 --repeat "$REPEAT" --csv "$X/mems/${c}_fb.csv"
      fi
    done
  done
  python3 "$S/summarize_mems.py" "$X"/mems/*.csv | tee "$X/summary.tsv"
  du -sh "$X"/*_t3 "$X"/*_t6
done
# All sizes in one table: individuals, then the columns of summarize_mems.py.
for size in ${SIZES//,/ }; do
  [ -s "$D/n$size/summary.tsv" ] || continue
  if [ ! -s "$D/summary_scaling.tsv" ] || [ "$size" = "${SIZES%%,*}" ]; then
    head -1 "$D/n$size/summary.tsv" | sed 's/^/individuals\t/' > "$D/summary_scaling.tsv"
  fi
  tail -n +2 "$D/n$size/summary.tsv" | sed "s/^/$size\t/" >> "$D/summary_scaling.tsv"
done
stamp "all done ($(free_gb "$DATA") GB free); results in $D/summary_scaling.tsv and $D/n*/"
