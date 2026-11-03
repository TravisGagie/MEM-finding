#!/bin/bash
# Zymo build, one capped step at a time (11 GB, no swap); stops at the first failure
cd ~/zymo/hpc; K=~/rz/MEM-finding/tool/bin/katka
step() { local name=$1; shift; echo "== $name $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$name: %e s wall, peak %M KB" "$@"
  local rc=$?; echo "$name exit $rc"; [ $rc = 0 ] || { echo "== STOPPED at $name"; exit 1; }; }
step prep $K prep --preset zymo ~/zymo/list.tsv zy; ls -l zy.S; df -h / | tail -1
step index-csa $K index --csa zy; ls -l zy.rlbwt zy.csa; df -h / | tail -1
echo "== done $(date +%T)"
