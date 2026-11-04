#!/bin/bash
# two-level locating (bb42c5a): index --twolevel-locate PREFIX, then find --twolevel --locate at L = 1000 and 100,
# pinned to CPU 9; checked against plain search (columns 1-5), ep - sp, and STPD over the parse with --locate
P=$1; W=$2; cd ~/tool/chr19; K=~/rz/MEM-finding/tool/bin/katka; mkdir -p $W
echo "== index --twolevel-locate $P $(date +%T)"
systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "tloc: %e s wall, peak %M KB" $K index --twolevel-locate $P > $W/tloc.$P.err 2>&1; echo "exit $?"; cat $W/tloc.$P.err; ls -l $P.tloc
for L in 1000 100; do
  echo "== L = $L $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "wall %e s, peak %M KB" taskset -c 9 $K find -L $L --format mems --csa --twolevel --locate --stats $P queries.fa $W/tlloc$L.$P.tsv > $W/tlloc$L.$P.err 2>&1; echo "exit $?"
  grep -E 'reads=|searches|locate|wall' $W/tlloc$L.$P.err
  python3 -I - $W/tlloc$L.$P.tsv $3/L$L.plain.tsv $3/L$L.stpd-parse-locate.tsv <<'PY'
import sys
t = [l.rstrip('\n').split('\t') for l in open(sys.argv[1])]
p = [l.rstrip('\n').split('\t') for l in open(sys.argv[2])]
s = [l.rstrip('\n').split('\t') for l in open(sys.argv[3])]
print('  %d MEMs (plain %d, STPD over the parse --locate %d)' % (len(t), len(p), len(s)))
print('  columns 1-5 equal plain search: %s' % (len(t) == len(p) and all(a[:5] == b[:5] for a, b in zip(t, p))))
print('  occurrences (column 6) equal ep - sp: %s' % all(int(a[5]) == int(a[4]) - int(a[3]) for a in t))
print('  occurrences and distinct sequences equal STPD over the parse --locate (per MEM): %s' % (len(t) == len(s) and all(a[:3] == b[:3] and a[5:7] == b[5:7] for a, b in zip(t, s))))
print('  occurrences in all %d, MEMs in more than one sequence %d' % (sum(int(a[5]) for a in t), sum(1 for a in t if int(a[6]) > 1)))
PY
done; echo "== done $(date +%T)"
