#!/bin/bash
# zymo preset: KeBaB's pipeline at the journal's fixed L = 25 for t = 10, 20, 30 (every MEM votes); and the check rule
# t = 1, T = 2 (runs e.2.1 / es.2.1 from 2026-11-03, rescored here with every MEM voting too)
cd ~/zymo/kebab2; K=~/rz/MEM-finding/tool/bin/katka; H=~/zymo/hpc/zy; R=~/zymo/all.fq; Z="-U -D 14 -c --csa -j 12 --stats"
KR=~/rz/MEM-finding/katka-readable/scripts; T=~/zymo/all.truth; B=~/zymo/hpc/zy.tbl
cap() { systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 "$@"; }
for t in 10 20 30; do
  echo "== f.$t $(date +%T)"
  cap /usr/bin/time -f "f.$t: %e s wall, %U s user, peak %M KB" $K find $Z -L 25 -K 20 --kebab-top $t --fixed-L $H $R f.$t > f.$t.err 2>&1; echo "exit $?"; grep -E 'kebab|wall,' f.$t.err
done
for r in f.10 f.20 f.30; do echo "-- $r: every MEM votes"; cap python3 -I $KR/credit_schemes.py $r $T $B | head -3; done
for r in e.2.1 es.2.1; do echo "-- $r: every MEM votes"; cap python3 -I $KR/credit_schemes.py <(zstdcat $r.zst) $T $B | head -3; done
for r in f.10 f.20 f.30; do zstd -q -T6 -19 --rm $r; done; echo "== done $(date +%T)"
