#!/bin/bash
# KeBaB Section 5.1 setting on Zymo's positive reads (no homopolymer compression, no low-complexity filter): KeBaB
# k = 20, one hash function, L = 40; ropebwt3 3.10 unmodified, one thread (kebab leaks every read with more than one thread; the paper also timed with one thread). (a) t longest pseudo-MEMs only;
# (b) T = 2t longest, whole-read fallback if fewer than t MEMs reach the (T+1)st pseudo-MEM; (c) whole reads.
cd ~/zymo/rb; B=bin; R=~/zymo/all.fq; S=~/research/mem-finding/tool/results/kebab_revision; L=40
cap() { systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 "$@"; }
tm() { local name=$1; shift; echo "== $name $(date +%T)"; /usr/bin/time -f "$name: %e s wall, %U s user, peak %M KB" "$@"; }
# resumed after the uncapped comparison (old longest.py, which held the whole-read MEMs in memory) was OOM-killed at
# 03:11; the t = 1 runs and (c) had finished and are reused; longest.py compare now streams in read order
echo "######## t = 1 (runs done before the stop; comparisons only)"
echo "-- (a) t = 1 against (c)"; cap python3 -I $S/longest.py compare a.1.mem c.mem.zst $R 1
echo "-- (b) t = 1 against (c)"; cap python3 -I $S/longest.py compare b.1.mem,b.1.fb.mem c.mem.zst $R 1
for f in a.1.fa a.1.mem b.1.all.fa b.1.fa b.1.mem b.1.fb.mem b.1.rerun.fq; do zstd -q -T6 --rm $f; done
for t in 5 10; do T=$((2 * t))
  echo "######## t = $t, T = $T"
  cap /usr/bin/time -f "a.$t scan: %e s wall, %U s user, peak %M KB" $B/kebab scan -i zymo.kbb -o a.$t.fa -l $L -s --top-t $t -t 1 $R 2> a.$t.scan.err; tail -1 a.$t.scan.err
  python3 -I $S/frag_stats.py $R < a.$t.fa
  cap /usr/bin/time -f "a.$t mem: %e s wall, %U s user, peak %M KB" $B/ropebwt3 mem -l $L -t 1 zymo.fmd a.$t.fa 2> a.$t.mem.err > a.$t.mem; tail -1 a.$t.mem.err
  cap /usr/bin/time -f "b.$t scan: %e s wall, %U s user, peak %M KB" $B/kebab scan -i zymo.kbb -o b.$t.all.fa -l $L -s --top-t $((T + 1)) -t 1 $R 2> b.$t.scan.err; tail -1 b.$t.scan.err
  /usr/bin/time -f "b.$t split: %e s wall" python3 -I $S/split_frags.py b.$t.all.fa $T b.$t.fa b.$t.lambda.tsv
  python3 -I $S/frag_stats.py $R < b.$t.fa
  cap /usr/bin/time -f "b.$t mem: %e s wall, %U s user, peak %M KB" $B/ropebwt3 mem -l $L -t 1 zymo.fmd b.$t.fa 2> b.$t.mem.err > b.$t.mem; tail -1 b.$t.mem.err
  /usr/bin/time -f "b.$t check: %e s wall" python3 -I $S/longest.py fallback b.$t.mem b.$t.lambda.tsv $t > b.$t.rerun.txt
  python3 -I - b.$t.rerun.txt $R > b.$t.rerun.fq <<'PY'
import sys
want = set(l.strip() for l in open(sys.argv[1]))
with open(sys.argv[2]) as F:
    while True:
        h = F.readline()
        if not h: break
        rec = h + F.readline() + F.readline() + F.readline()
        if h[1:].split()[0] in want: sys.stdout.write(rec)
PY
  cap /usr/bin/time -f "b.$t fallback mem: %e s wall, %U s user, peak %M KB" $B/ropebwt3 mem -l $L -t 1 zymo.fmd b.$t.rerun.fq 2> b.$t.fb.err > b.$t.fb.mem; tail -1 b.$t.fb.err
  echo "-- (a) t = $t against (c)"; cap python3 -I $S/longest.py compare a.$t.mem c.mem.zst $R $t
  echo "-- (b) t = $t against (c)"; cap python3 -I $S/longest.py compare b.$t.mem,b.$t.fb.mem c.mem.zst $R $t
  for f in a.$t.fa a.$t.mem b.$t.all.fa b.$t.fa b.$t.mem b.$t.fb.mem b.$t.rerun.fq; do zstd -q -T6 --rm $f; done
done
ls -l; df -h / | tail -1; echo "== done $(date +%T)"
