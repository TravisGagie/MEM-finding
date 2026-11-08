#!/bin/bash
cd ~/zymo/newrun/reg; S=~/research/mem-finding/katka-readable/scripts; T=~/zymo/all.truth; B=~/zymo/hpc/zy.tbl
echo "== L30 credit_schemes $(date +%T)"; python3 -I $S/credit_schemes.py <(zstdcat new.L30.zst) $T $B
echo "== L30 score_spu2 $(date +%T)";     python3 -I ~/zymo/score_spu2.py <(zstdcat new.L30.zst) $T $B
echo "== X credit_schemes $(date +%T)";   python3 -I $S/credit_schemes.py <(zstdcat new.X.zst) $T $B
echo "== X score_spu2 $(date +%T)";       python3 -I ~/zymo/score_spu2.py <(zstdcat new.X.zst) $T $B
echo "== done $(date +%T)"
