#!/bin/bash
# c33: experiments/stpd_parse_chr19.sh (STPD over the parse, with and without --locate; character-level STPD through
# the parse oracle; two-level with blocks; plain search), plus katka's STPD over the text (PREFIX.S in memory) with
# the same binary and sample, in the same way (pinned to CPU 9, 11 GB cgroup, MEMs compared with plain search).
set -uo pipefail
cd ~/research/rlbwt-tricks/chr19; K=~/research/mem-finding/tool/bin/katka; W=$HOME/research/rlbwt-tricks/chr19/stpd_parse33; PAT=$HOME/research/rlbwt-tricks/chr19/queries.fa
[ -s c33.stpd ] && mv c33.stpd c33.stpd.old          # rebuilt by the script with the current binary
rm -f c33.ptxt c33.pstpd c33.rptxt
mkdir -p $W; echo "== index --parse-oracle c33 (the dictionary at 2 bits per base) $(date +%T)"
systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "parse-oracle: %e s wall, peak %M KB" $K index --parse-oracle c33 > $W/build_ptxt.log 2>&1; echo "exit $?"
grep -iE "parse oracle|dictionary|exception|wall" $W/build_ptxt.log; ls -l c33.ptxt
PIN=9 ~/research/mem-finding/tool/experiments/stpd_parse_chr19.sh $HOME/research/rlbwt-tricks/chr19/c33 $PAT $W; echo "script exit $?"
SUM=$W/summary.txt
for L in 1000 100; do
  for v in "stpd-text --stpd --no-intervals" "stpd-text-locate --stpd --locate"; do set -- $v; name=L$L.$1; shift
    echo "[$(date +%T)] $name"
    systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "wall %e s, peak memory %M KB" taskset -c 9 $K find --format mems --stats -L $L "$@" c33 $PAT $W/$name.tsv > $W/$name.log 2>&1
    echo "  $(printf '%-26s' $name) $(grep -o ' MEMs=[0-9]*' $W/$name.log | head -1 | tr -d ' '), search $(grep -o 'wall [0-9.]* s' $W/$name.log | head -1); $(tail -1 $W/$name.log)" | tee -a $SUM
    if cmp -s <(cut -f1-3 $W/$name.tsv) <(cut -f1-3 $W/L$L.plain.tsv); then echo "      $name: the same MEMs as L$L.plain" | tee -a $SUM; else echo "      $name: MEMS DIFFER from L$L.plain" | tee -a $SUM; fi
  done
done
ls -l c33.* | tee -a $SUM; echo "== done $(date +%T)"
