#!/bin/bash
cd ~/zymo/newrun/reg; echo "== maxmem_roc (streamed), longest MEMs with KeBaB $(date +%T)"
python3 -I ~/research/mem-finding/tool/scripts/maxmem_roc_stream.py <(zstdcat new.X.zst) ~/zymo/all.truth ~/zymo/hpc/zy.tbl ~/zymo/all.fq
echo "== done $(date +%T)"
