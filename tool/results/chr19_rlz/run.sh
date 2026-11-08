#!/bin/bash
cd ~/research/rlbwt-tricks/chr19/rlz; K=~/research/mem-finding/tool/bin/katka
echo "== index --parse-oracle $(date +%T)"
/usr/bin/time -f "parse-oracle: %e s wall, %U s user, peak %M KB" $K index --parse-oracle c33 2>&1 | tail -3
ls -l c33.ptxt
echo "== rlz_parse $(date +%T)"
/usr/bin/time -f "rlz_parse: %e s wall, %U s user, peak %M KB" ~/research/mem-finding/tool/bin/rlz_parse c33.ptxt 61707364 2>&1
echo "== done $(date +%T)"
