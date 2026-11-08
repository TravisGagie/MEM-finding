#!/bin/bash
# c33 with a sparser parse: closed syncmers k = 200, s = 20 (the text, BWT and character-level STPD sample are c33's,
# by symbolic link; the new prefix is c). Then experiments/stpd_parse_chr19.sh as for c33 at 51,12.
cd ~/research/rlbwt-tricks/chr19/sparse; K=~/research/mem-finding/tool/bin/katka
echo "== index --preset chr19 --syncmer 200,20 $(date +%T)"
systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "index: %e s wall, peak %M KB" $K index --preset chr19 --syncmer 200,20 c > index.err 2>&1; echo "exit $?"; cat index.err; ls -l c.*
PIN=9 ~/research/mem-finding/tool/experiments/stpd_parse_chr19.sh $HOME/research/rlbwt-tricks/chr19/sparse/c $HOME/research/rlbwt-tricks/chr19/queries.fa $HOME/research/rlbwt-tricks/chr19/sparse/work; echo "script exit $?"
ls -l c.*; df -h / | tail -1; echo "== done $(date +%T)"
