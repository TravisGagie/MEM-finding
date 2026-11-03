#!/bin/bash
cd ~/tool/chr19; K=~/rz/MEM-finding/tool/bin/katka; X=~/rz/MEM-finding/tool/experiments/stpd_index_chr19.sh
if [ ! -s c33.stpd ]; then echo "== index --stpd c33 $(date +%T)"
  /usr/bin/time -f "index-stpd: %e s wall, %U s user, peak %M KB" $K index --stpd c33 > index-stpd.33.err 2>&1; rc=$?; tail -4 index-stpd.33.err; echo "exit $rc"; [ $rc = 0 ] || exit 1
fi; ls -l c33.*
for L in 1000 100; do
  echo "== comparison L = $L $(date +%T)"
  BIGBWT=~/rz/MEM-finding/katka/deps/Big-BWT/bigbwt OURS=$PWD/c33 PIN=9 $X ~/rz/STPD-index $PWD/list33abs.tsv $PWD/queries.fa $PWD/stpd33 $L > stpd33.L$L.log 2>&1; rc=$?
  echo "exit $rc"; cp stpd33/summary.txt stpd33/summary.L$L.txt; cat stpd33/summary.L$L.txt; df -h / | tail -1; [ $rc = 0 ] || exit 1
done; echo "== done $(date +%T)"
