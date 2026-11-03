#!/bin/bash
# KeBaB stopping rules on Zymo (tool/NOTES.md, "KeBaB stopping rules for the journal revision"), tool at 7500ed0
cd ~/zymo/kebab; K=~/rz/MEM-finding/tool/bin/katka; H=~/zymo/hpc/zy; R=~/zymo/all.fq; Z="-U -D 14 -c --csa -j 12 --stats"
S=~/rz/MEM-finding/tool/scripts; T=~/zymo/all.truth; B=~/zymo/hpc/zy.tbl
step() { local name=$1; shift; echo "== $name $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$name: %e s wall, %U s user, peak %M KB" "$@" > $name.err 2>&1
  local rc=$?; grep -E 'reads=|kebab|searches|wall,' $name.err; echo "$name exit $rc"; [ $rc = 0 ] || { echo "== STOPPED at $name"; exit 1; }; }
for t in 1 3 5; do
  step a.$t  $K find $Z -L 20 -K 20 --kebab-top $t   --read-stats a.$t.tsv  $H $R a.$t
  step b.$t  $K find $Z -L 20 -K 20 --kebab-batch $t --read-stats b.$t.tsv  $H $R b.$t
  step b1.$t $K find $Z -L 1  -K 20 --kebab-batch $t --read-stats b1.$t.tsv $H $R b1.$t
  step c.$t  $K find $Z -L 20 -K 20 --longest $t     --read-stats c.$t.tsv  $H $R c.$t
  step d.$t  $K find $Z -L 1        --longest $t                            $H $R d.$t
done
echo "== accuracy (maxmem_roc, longest MEMs vote) $(date +%T)"
for t in 1 3 5; do for r in a b b1 c d; do
  echo "-- $r.$t"; python3 -I $S/maxmem_roc_stream.py $r.$t $T $B $R | sed -n '3p'; done; done
echo "== comparisons $(date +%T)"
for t in 1 3 5; do
  for p in "a c" "b c" "a d" "b1 d" "c d"; do set -- $p
    echo "-- $1.$t against $2.$t"; python3 -I $S/compare_longest_t.py $1.$t $2.$t $R $t; done; done
echo "== compressing $(date +%T)"; for f in [abcd].[135] b1.[135]; do zstd -q -T6 -19 --rm $f; done; ls -l
echo "== done $(date +%T)"
