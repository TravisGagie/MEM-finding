#!/bin/bash
# KeBaB Section 5.1 setting on Zymo's positive reads (no homopolymer compression, no low-complexity filter): KeBaB
# k = 20, one hash function, L = 40; ropebwt3 3.10 unmodified, one thread (kebab leaks every read with more than one thread; the paper also timed with one thread). (a) t longest pseudo-MEMs only;
# (b) T = 2t longest, whole-read fallback if fewer than t MEMs reach the (T+1)st pseudo-MEM; (c) whole reads.
cd ~/zymo/rb; B=bin; R=~/zymo/all.fq; S=~/rz/MEM-finding/tool/results/kebab_revision; L=40
cap() { systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 "$@"; }
tm() { local name=$1; shift; echo "== $name $(date +%T)"; /usr/bin/time -f "$name: %e s wall, %U s user, peak %M KB" "$@"; }
# the Section 5.1 check rule (b) with T = t + 1 and T = t + 2 (t = 1, T = 2 is the T = 2t run of run.sh)
for tT in "10 15"; do set -- $tT; t=$1; T=$2
  echo "######## t = $t, T = $T (check rule only)"
  cap /usr/bin/time -f "b.$t.T$T scan: %e s wall, %U s user, peak %M KB" $B/kebab scan -i zymo.kbb -o b.$t.T$T.all.fa -l $L -s --top-t $((T + 1)) -t 1 $R 2> b.$t.T$T.scan.err; tail -1 b.$t.T$T.scan.err
  /usr/bin/time -f "b.$t.T$T split: %e s wall" python3 -I $S/split_frags.py b.$t.T$T.all.fa $T b.$t.T$T.fa b.$t.T$T.lambda.tsv
  python3 -I $S/frag_stats.py $R < b.$t.T$T.fa
  cap /usr/bin/time -f "b.$t.T$T mem: %e s wall, %U s user, peak %M KB" $B/ropebwt3 mem -l $L -t 1 zymo.fmd b.$t.T$T.fa 2> b.$t.T$T.mem.err > b.$t.T$T.mem; tail -1 b.$t.T$T.mem.err
  /usr/bin/time -f "b.$t.T$T check: %e s wall" python3 -I $S/longest.py fallback b.$t.T$T.mem b.$t.T$T.lambda.tsv $t > b.$t.T$T.rerun.txt
  python3 -I - b.$t.T$T.rerun.txt $R > b.$t.T$T.rerun.fq <<'PY'
import sys
want = set(l.strip() for l in open(sys.argv[1]))
with open(sys.argv[2]) as F:
    while True:
        h = F.readline()
        if not h: break
        rec = h + F.readline() + F.readline() + F.readline()
        if h[1:].split()[0] in want: sys.stdout.write(rec)
PY
  cap /usr/bin/time -f "b.$t.T$T fallback mem: %e s wall, %U s user, peak %M KB" $B/ropebwt3 mem -l $L -t 1 zymo.fmd b.$t.T$T.rerun.fq 2> b.$t.T$T.fb.err > b.$t.T$T.fb.mem; tail -1 b.$t.T$T.fb.err
  echo "-- (b) t = $t against (c)"; cap python3 -I $S/longest.py compare b.$t.T$T.mem,b.$t.T$T.fb.mem c.mem.zst $R $t
  for f in b.$t.T$T.all.fa b.$t.T$T.fa b.$t.T$T.mem b.$t.T$T.fb.mem b.$t.T$T.rerun.fq; do zstd -q -T6 --rm $f; done
done
ls -l; df -h / | tail -1; echo "== done $(date +%T)"
