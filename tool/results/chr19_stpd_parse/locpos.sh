#!/bin/bash
# position-by-position comparison (65e46c4): two-level locating against STPD over the parse, --locate-out, cmp
cd ~/tool/chr19; K=~/rz/MEM-finding/tool/bin/katka; W=~/tool/chr19/locpos; mkdir -p $W
for P in c33 c65; do n=${P#c}
  for L in 1000 100; do
    echo "== $P L = $L $(date +%T)"
    systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "two-level --locate: wall %e s, peak %M KB" taskset -c 9 $K find -L $L --format mems --csa --twolevel --locate --locate-out $W/tlpos$L.$n.txt $P queries.fa $W/tlloc$L.$n.tsv > $W/tlloc$L.$n.err 2>&1; echo "exit $?"
    systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "STPD over the parse --locate: wall %e s, peak %M KB" taskset -c 9 $K find -L $L --format mems --csa --stpd-parse --locate --locate-out $W/sppos$L.$n.txt $P queries.fa $W/sploc$L.$n.tsv > $W/sploc$L.$n.err 2>&1; echo "exit $?"
    grep -hE 'reads=|wall,' $W/tlloc$L.$n.err $W/sploc$L.$n.err
    if cmp $W/tlpos$L.$n.txt $W/sppos$L.$n.txt; then echo "  $P L = $L: positions byte-identical ($(wc -l < $W/tlpos$L.$n.txt) MEMs, $(stat -c %s $W/tlpos$L.$n.txt) bytes)"; else echo "  $P L = $L: POSITIONS DIFFER"; fi
    head -c 300 $W/tlpos$L.$n.txt | head -2 | cut -c1-150
  done
done; df -h / | tail -1; echo "== done $(date +%T)"
