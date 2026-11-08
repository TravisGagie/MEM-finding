#!/bin/bash
cd ~/zymo/newrun/reg; echo "== same longest MEMs, KeBaB vs no KeBaB $(date +%T)"
python3 -I ~/research/mem-finding/tool/scripts/same_longest_mems.py <(zstdcat new.X.zst) <(zstdcat new.X.nokebab.zst) ~/zymo/all.fq
echo "== maxmem_roc on the no-KeBaB run $(date +%T)"
python3 -I ~/research/mem-finding/tool/scripts/maxmem_roc_stream.py <(zstdcat new.X.nokebab.zst) ~/zymo/all.truth ~/zymo/hpc/zy.tbl ~/zymo/all.fq | head -4
echo "== done $(date +%T)"
