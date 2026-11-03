#!/bin/bash
cd ~/tool/silva/blocks; K=~/rz/MEM-finding/tool/bin/katka; B=~/rz/MEM-finding/tool/bin/katka-blocks; P=~/tool/silva/bac
systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 $K find -j 1 --preset silva $P reads_u.fq listing_u.out 2> listing_u.err; tail -1 listing_u.err
for ks in 8,12 16; do echo "== k = $ks $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "k=$ks: %e s wall, %U s user, peak %M KB" $B $P reads_u.fq listing_u.out -k $ks -L 30 -n 8000 > blocks_u_k$ks.txt 2>&1; echo "exit $?"; done
echo "== done $(date +%T)"
