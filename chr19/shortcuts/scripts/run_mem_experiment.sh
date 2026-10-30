#!/bin/bash
# Times MEM finding with two indexes (query_mems), with and without
# shortcuts, for several configurations.
#
#   run_mem_experiment.sh INDEX REVERSED_INDEX PATTERN_DIR OUTDIR
#
# INDEX             directory with rlbwt.bin and ell<L>_theta<T>.shortcuts
#                   for T (build_index)
# REVERSED_INDEX    the same for the reversed collection (build_index --reversed)
# PATTERN_DIR       directory with patterns_<p>.txt (generate_patterns.py)
#
# Environment variables:
#   CONFIGS=ell16_theta0,ell25_theta0,ell50_theta0,ell100_theta0
#   PERCENTS=0,0.1,1   REPEAT=3
#   RUNS="fb:1 bml:50,100"   algorithm:minimum lengths, one query_mems call each
#   PIN=<cpu>        run the queries on this CPU only (with taskset)
#   BIN=<directory with query_mems>
#
# Writes OUTDIR/mems/<config>_<p>_<algorithm>.csv, OUTDIR/mem_results.csv
# and OUTDIR/mem_summary.txt.
set -euo pipefail

if [ $# -ne 4 ]; then
  sed -n '2,22p' "$0"
  exit 1
fi
INDEX=$1
REVERSED=$2
PATTERNS=$3
OUT=$4
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-$HERE/../build}
CONFIGS=${CONFIGS:-ell16_theta0,ell25_theta0,ell50_theta0,ell100_theta0}
PERCENTS=${PERCENTS:-0,0.1,1}
REPEAT=${REPEAT:-3}
RUNS=${RUNS:-"fb:1 bml:50,100"}
RUN=()
if [ -n "${PIN:-}" ]; then RUN=(taskset -c "$PIN"); fi
stamp() { echo "[$(date '+%F %T')] $*"; }

mkdir -p "$OUT/mems"
for config in ${CONFIGS//,/ }; do
  for p in ${PERCENTS//,/ }; do
    for run in $RUNS; do
      algorithm=${run%%:*}
      lengths=${run#*:}
      stamp "$config, $p%, $algorithm"
      ${RUN[@]+"${RUN[@]}"} "$BIN/query_mems" --index "$INDEX" --reversed-index "$REVERSED" --shortcuts "$config" \
        --patterns "$PATTERNS/patterns_$p.txt" --algorithm "$algorithm" --min-length "$lengths" \
        --repeat "$REPEAT" --csv "$OUT/mems/${config}_${p}_${algorithm}.csv"
    done
  done
done

stamp "collecting results"
python3 "$HERE/analyze_mems.py" "$OUT" "$INDEX" "$REVERSED"
stamp "done"
