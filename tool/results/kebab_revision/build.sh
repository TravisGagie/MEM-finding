#!/bin/bash
# KeBaB Section 5.1 setting: a ropebwt3 index of the uncompressed Zymo reference (both strands) and KeBaB's filter
# (k = 20, eps = 0.1, one hash function, as in the paper's appendix). Each step capped at 11 GB, no swap.
cd ~/zymo/rb; X="/run/media/travis/Travis backup/zymo_plain"
step() { local name=$1; shift; echo "== $name $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$name: %e s wall, %U s user, peak %M KB" "$@" > $name.err 2>&1
  local rc=$?; tail -3 $name.err; echo "$name exit $rc"; [ $rc = 0 ] || { echo "== STOPPED at $name"; exit 1; }; }
step ropebwt3-build bin/ropebwt3 build -t 12 -m 1g -do zymo.fmd $(cut -f2 ~/zymo/list.tsv)
ls -l zymo.fmd
echo "== uncompressed reference for kebab build (external drive) $(date +%T)"
for f in $(cut -f2 ~/zymo/list.tsv); do zcat $f; done > "$X/ref.fa"; ls -l "$X/ref.fa"
step kebab-build bin/kebab build -k 20 -e 0.1 -f 1 -t 12 -o zymo "$X/ref.fa"
ls -l zymo.*; df -h / | tail -1; echo "== done $(date +%T)"
