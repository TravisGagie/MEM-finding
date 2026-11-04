#!/bin/bash
# Cliffy's approximate document listing, scored on our reads with one command (tool/NOTES.md, "Cliffy").
#
#   tool/scripts/score_cliffy_listing.sh [RESULTS_DIR]
#
# RESULTS_DIR is the directory katka/silva/run_cliffy.sh wrote (default: the newest
# katka/silva/work/compare/cliffy-*), whose runs/ holds Cliffy's listings: runs/<index>.<region>_<mate>.p1.all,
# from `cliffy run --taxcomp --ftab` (index full_text or minimizers).  For each index and region it writes
# RESULTS_DIR/scores/<index>.<region>.listing.csv with tool/scripts/cliffy_listing_score.py (their "method2",
# the scoring of approximate listing in Cliffy's paper) and prints a table of genus-level accuracies.
# Needs the files run_cliffy.sh used: katka/silva/work/ref, katka/silva/work/reads, the SILVA taxonomy and
# katka/cliffy-experiments (their code, for loading the files).  Environment: S (default katka/silva),
# EXP (cliffy-experiments/src), TAX, REGIONS (default "V1_V2 V3_V4 V4_V4 V4_V5"), INDEXES (default
# "full_text minimizers"), THREADS.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); TOP=$(cd "$HERE/../.." && pwd)
S=${S:-$TOP/katka/silva}; W=$S/work
EXP=${EXP:-$(dirname "$S")/cliffy-experiments/src}; TAX=${TAX:-$S/exp1_data/tax_slv_ssu_138.1.txt}
REGIONS=${REGIONS:-"V1_V2 V3_V4 V4_V4 V4_V5"}; INDEXES=${INDEXES:-"full_text minimizers"}; THREADS=${THREADS:-$(nproc)}
R=${1:-$(ls -dt $W/compare/cliffy-* 2> /dev/null | head -1)}
[ -n "$R" ] && [ -d "$R/runs" ] || { echo "no Cliffy results (run katka/silva/run_cliffy.sh, or give RESULTS_DIR)"; exit 1; }
[ -f "$EXP/classify_silva_readset.py" ] || { echo "no cliffy-experiments at $EXP (git clone https://github.com/oma219/cliffy-experiments; set EXP)"; exit 1; }
for f in "$TAX" "$W/ref/rna/doc_to_traversal.txt" "$W/ref/trav_to_length.txt"; do [ -s "$f" ] || { echo "missing $f"; exit 1; }; done
mkdir -p "$R/scores"
printf "%-12s" index; for r in $REGIONS; do printf " %8s" $r; done; printf " %8s\n" mean
for c in $INDEXES; do
  line=$(printf "%-12s" $c); sum=0; k=0
  for r in $REGIONS; do
    m1=$R/runs/$c.${r}_1.p1.all m2=$R/runs/$c.${r}_2.p1.all o=$R/scores/$c.$r.listing.csv
    if [ ! -s "$m1" ] || [ ! -s "$m2" ]; then line+=$(printf " %8s" "-"); continue; fi
    [ -s "$o" ] || python3 "$HERE/cliffy_listing_score.py" --mate1-listings "$m1" --mate2-listings "$m2" \
        --doc-id-to-traversal $W/ref/rna/doc_to_traversal.txt --silva-tax-ranks "$TAX" \
        --readset-truthset $W/reads/aquatic/$r/${r}_seqtax.txt --trav-to-length $W/ref/trav_to_length.txt \
        --exp-src "$EXP" --output "$o" --threads $THREADS > "$o.log" 2>&1 || { echo "scoring $c $r failed (see $o.log)"; exit 1; }
    a=$(grep "^cliffy_listing,genus," "$o" | grep -o "theirs=[0-9.]*" | cut -d= -f2)
    line+=$(printf " %7.2f%%" $(echo "100 * $a" | bc -l)); sum=$(echo "$sum + $a" | bc -l); k=$((k + 1))
  done
  [ $k -gt 0 ] && line+=$(printf " %7.2f%%" $(echo "100 * $sum / $k" | bc -l))
  echo "$line"
done
echo "genus level, Cliffy's rule (theirs= in each csv; the files also give the other levels, and accuracy counting"
echo "reads without a match as wrong, and strict, counting ties as wrong).  Scores: $R/scores/*.listing.csv"
