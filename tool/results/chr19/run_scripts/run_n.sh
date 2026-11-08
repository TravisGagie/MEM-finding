#!/bin/bash
# chr19 with the new tool (tool/NOTES.md, "Chromosome 19"), one collection size: $1 = 33, 65 or 129 sequences
set -o pipefail
n=$1; H=$HOME/research/rlbwt-tricks/shortcuts/chr19_scaling/hprc_chr19; cd ~/research/rlbwt-tricks/chr19; K=~/research/mem-finding/tool/bin/katka; BB=~/research/mem-finding/katka/deps/Big-BWT/bigbwt
step() { local name=$1; shift; echo "== $name $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$name: %e s wall, %U s user, peak %M KB" "$@" > $name.$n.err 2>&1
  local rc=$?; grep -vE '^\s*$' $name.$n.err | tail -4; echo "$name exit $rc"; [ $rc = 0 ] || { echo "== STOPPED at $name"; exit 1; }; }
if [ ! -s c$n.rlbwt ]; then
  mkdir -p fa$n; head -n $n $H/collection.txt | awk -v d=fa$n '{f = sprintf("%s/s%03d.fa", d, NR); print ">s" NR > f; print > f; close(f)}'
  ls fa$n/*.fa | sed 's/^/-\t/' > list$n.tsv
  step prep $K prep list$n.tsv c$n; rm -rf fa$n; ls -l c$n.S
  step bigbwt $BB c$n.S; mv c$n.S.bwt c$n.bwt; rm -f c$n.S.log c$n.S.bigbwt.log
  step index $K index --preset chr19 c$n; ls -l c$n.*; df -h / | tail -1
fi
[ -s queries.fa ] || tr -c 'ACGT\n' '\n' < $H/queries.txt | awk 'length($0) > 0 { print ">p" NR - 1; print }' > queries.fa
for L in 1000 100; do
  step bml$L taskset -c 9 $K find --preset chr19 -L $L --stats c$n queries.fa bml$L.$n.tsv
  step nj$L taskset -c 9 $K find --preset chr19 -L $L --no-jumps --stats c$n queries.fa nj$L.$n.tsv
  step plain$L taskset -c 9 $K find -L $L --format mems --csa --stats c$n queries.fa plain$L.$n.tsv
  cmp bml$L.$n.tsv plain$L.$n.tsv && cmp nj$L.$n.tsv plain$L.$n.tsv && echo "L = $L: two-level (blocks, no jumps) = plain, $(wc -l < plain$L.$n.tsv) MEMs"
done
step fb1000 taskset -c 9 $K find --preset chr19 --fb --stats c$n queries.fa fb1000.$n.tsv
cmp fb1000.$n.tsv plain1000.$n.tsv && echo "fb = plain (L = 1000)"
echo "== done $(date +%T)"
