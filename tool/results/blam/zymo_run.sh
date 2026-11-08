#!/bin/bash
# B_lambda on Zymo (3ab8589): per-species tags, lambda = 0, 10, 30, 100, 1000; text and BWT and temporary files on the
# external drive; capped at 11 GB. The run count and the tag runs it writes must equal the index's zy.sp.tagruns.
X="/run/media/travis/Travis backup"; P="$X/zymo_bwt/zy"; K=~/research/mem-finding/tool/bin/katka
echo "== B_lambda $(date +%T)"
systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "blam: %e s wall, %U s user, peak %M KB" \
  $K index --tags --by-tag --no-rmq --tag-out /dev/null --tag-lcp 0,10,30,100,1000 --tmp "$X/blam_tmp" "$P" > ~/zymo/blam.zymo.err 2>&1; echo "exit $?"
grep -vE '^\s*$' ~/zymo/blam.zymo.err | tail -25
ls -l "$X/zymo_bwt/"
cmp "$P.sp.tagruns" ~/zymo/hpc/zy.sp.tagruns && echo "tag runs: byte-identical to the index's zy.sp.tagruns"
echo "== done $(date +%T)"
