#!/bin/bash
# convert KATKA's tag runs, then check them against the BWT; capped at 11 GB, stops at the first failure
cd ~/zymo/hpc
step() { local name=$1; shift; echo "== $name $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$name: %e s wall, peak %M KB" "$@"
  local rc=$?; echo "$name exit $rc"; [ $rc = 0 ] || { echo "== STOPPED at $name"; exit 1; }; }
step convert python3 -I ~/zycheck/scripts/convert_tagruns.py zy.sp.tagruns.katka zy.sp.tagruns.conv; ls -l zy.sp.tagruns.conv; df -h / | tail -1
step check ~/zycheck/conv/check_tagruns zy zy.sp.tagruns.conv 5000000 1
echo "== done $(date +%T)"
