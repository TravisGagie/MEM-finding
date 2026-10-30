#!/bin/bash
# Runs the experiments of the paper on one collection.
#
#   run_experiment.sh COLLECTION REFERENCE_LENGTH OUTDIR [BWT]
#
# COLLECTION        one sequence per line (see prepare_collection.py)
# REFERENCE_LENGTH  R is the first REFERENCE_LENGTH characters of the text
# BWT               optional: a BWT of COLLECTION from Big-BWT, for
#                   collections too large for libsais in memory
#
# Environment variables (defaults are those of the paper's cere runs):
#   ELLS=25,50,100,200   THETAS=0,100,...,1000   PERCENTS=0,1   REPEAT=5
#   BIN=<directory with build_index and query_index>
#   PIN=<cpu>        run the queries on this CPU only (with taskset), so
#                    that all timings come from the same kind of core
#   QUERY_ONLY=1     reuse OUTDIR/index and OUTDIR/patterns_*.txt
#   RLBWT=<file>     reuse an RLBWT built earlier (build_index --rlbwt)
#
# Writes OUTDIR/index/ (the RLBWT and shortcuts, with build.csv),
# OUTDIR/patterns_*.txt, OUTDIR/queries/*.csv, OUTDIR/results.csv and
# OUTDIR/summary.txt.
set -euo pipefail

if [ $# -lt 3 ]; then
  sed -n '2,21p' "$0"
  exit 1
fi
COLLECTION=$1
REFERENCE_LENGTH=$2
OUT=$3
BWT=${4:-}
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-$HERE/../build}
ELLS=${ELLS:-25,50,100,200}
THETAS=${THETAS:-0,100,200,300,400,500,600,700,800,900,1000}
PERCENTS=${PERCENTS:-0,1}
REPEAT=${REPEAT:-5}

mkdir -p "$OUT/queries"
stamp() { echo "[$(date '+%F %T')] $*"; }
RUN=()
if [ -n "${PIN:-}" ]; then RUN=(taskset -c "$PIN"); fi

if [ -z "${QUERY_ONLY:-}" ]; then
stamp "building the index"
rm -f "$OUT/index/build.csv"
if [ -n "${RLBWT:-}" ]; then
  "$BIN/build_index" --text "$COLLECTION" --rlbwt "$RLBWT" --reference-length "$REFERENCE_LENGTH" \
    --out "$OUT/index" --ell "$ELLS" --theta "$THETAS"
elif [ -n "$BWT" ]; then
  "$BIN/build_index" --text "$COLLECTION" --bwt "$BWT" --reference-length "$REFERENCE_LENGTH" \
    --out "$OUT/index" --ell "$ELLS" --theta "$THETAS"
else
  "$BIN/build_index" --text "$COLLECTION" --reference-length "$REFERENCE_LENGTH" \
    --out "$OUT/index" --ell "$ELLS" --theta "$THETAS"
fi

stamp "generating patterns"
python3 "$HERE/generate_patterns.py" "$COLLECTION" "$OUT" "$PERCENTS"
fi

stamp "querying"
for shortcuts in "$OUT"/index/*.shortcuts; do
  config=$(basename "$shortcuts" .shortcuts)
  for p in ${PERCENTS//,/ }; do
    ${RUN[@]+"${RUN[@]}"} "$BIN/query_index" --rlbwt "$OUT/index/rlbwt.bin" --shortcuts "$shortcuts" \
      --patterns "$OUT/patterns_$p.txt" --repeat "$REPEAT" --csv "$OUT/queries/${config}_$p.csv"
  done
done

stamp "collecting results"
python3 "$HERE/collect_results.py" "$OUT"
python3 "$HERE/analyze.py" "$OUT"
stamp "done"
