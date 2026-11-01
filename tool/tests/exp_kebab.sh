#!/bin/bash
# The experiment of Section 3.2 ("KeBaB inside BML"): on queries mostly absent from the text, with a fixed,
# large L, compare
#   KeBaB then BML      KeBaB's pseudo-MEMs (every k-mer of the query hashed), then BML on those of length >= L
#   plain BML
#   KeBaB inside BML    each window's k-mers tested from its right end, stopping at the first that fails
# single-threaded, on the same index (RLCSA), with the median of R runs; checks that all three give the same
# MEMs, and reports k-mer tests and backward steps.
#   exp_kebab.sh PREFIX queries.fa [L=1000] [k=20] [R=5]
# PREFIX needs PREFIX.csa and PREFIX.k<k>.kbf (katka index --csa --kebab k).  tests/gen_novel.py makes
# synthetic queries: random sequence with a few segments copied from the collection.
set -euo pipefail
P=$1; Q=$2; L=${3:-1000}; k=${4:-20}; R=${5:-5}
K=$(cd "$(dirname "$0")/.." && pwd)/bin/katka
W=$(mktemp -d)
run() {   # name, options
  local name=$1; shift
  local times=()
  for ((i = 0; i < R; ++i)); do
    "$K" find -L "$L" --format mems --csa --stats "$@" "$P" "$Q" "$W/$name.tsv" 2> "$W/$name.log"
    times+=("$(grep -o 'wall [0-9.]*' "$W/$name.log" | cut -d' ' -f2)")
  done
  local med; med=$(printf '%s\n' "${times[@]}" | sort -g | awk '{a[NR]=$1} END {print (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2}')
  local st; st=$(grep searches "$W/$name.log")
  printf '%-18s %8.3f s   k-mer tests %10s   restarts %8s   backward steps %10s   MEMs %s\n' "$name" "$med" \
    "$(echo "$st" | grep -o 'kmer_tests=[0-9]*' | cut -d= -f2)" "$(echo "$st" | grep -o 'kmer_restarts=[0-9]*' | cut -d= -f2)" \
    "$(echo "$st" | grep -o 'char_steps=[0-9]*' | cut -d= -f2)" "$(wc -l < "$W/$name.tsv")"
}
echo "L = $L, k = $k, $R runs each (median), queries: $(grep -c '>' "$Q") sequences, $(grep -v '>' "$Q" | tr -d '\n' | wc -c) characters"
run kebab-then-bml -K "$k"
run plain-bml
run kebab-inside-bml --inline-kebab "$k"
cmp -s "$W/plain-bml.tsv" "$W/kebab-then-bml.tsv" && cmp -s "$W/plain-bml.tsv" "$W/kebab-inside-bml.tsv" && echo "all three: the same MEMs" || { echo "DIFFERENT MEMs"; exit 1; }
rm -rf "$W"
