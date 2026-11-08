#!/bin/bash
cd ~/zymo/hpc; K=~/research/mem-finding/tool/bin/katka
echo "== index-kebab $(date +%T)"
/usr/bin/time -f "index-kebab: %e s wall, peak %M KB" $K index --kebab 20 zy; rc=$?; echo "index-kebab exit $rc"
[ $rc = 0 ] || { echo "== STOPPED at index-kebab"; exit 1; }
ls -l zy.k20.kbf katka/zy.k20.kbf
cmp zy.k20.kbf katka/zy.k20.kbf && echo "KeBaB filter: byte-identical to KATKA's" || echo "KeBaB filter: DIFFERENT from KATKA's"
echo "== done $(date +%T)"
