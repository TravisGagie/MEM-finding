#!/bin/bash
cd ~/zymo/hpc; K=~/rz/MEM-finding/tool/bin/katka
step() { local name=$1; shift; echo "== $name $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$name: %e s wall, peak %M KB" "$@"
  local rc=$?; echo "$name exit $rc"; [ $rc = 0 ] || { echo "== STOPPED at $name"; exit 1; }; }
step index-tags $K index --tags --by-tag --no-rmq zy; ls -l zy.tag zy.sp.tagruns; df -h / | tail -1
step info $K info zy
echo "== done $(date +%T)"
