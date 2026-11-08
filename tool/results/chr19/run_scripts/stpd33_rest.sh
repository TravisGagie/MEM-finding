#!/bin/bash
# The rest of experiments/stpd_index_chr19.sh (c09732c) for 16 individuals, without their forward-backward run
# (--fb), which ran for 54 min without finishing on L = 1000 (their BML took 45 s): the same commands for the other
# runs, and the script's comparison code, extracted verbatim (stpd33_compare.part). Their index is reused.
set -uo pipefail
W=$HOME/research/rlbwt-tricks/chr19/stpd33; cd $W; K=$HOME/research/mem-finding/tool/bin/katka; DRV=$HOME/research/mem-finding/tool/bin/stpd-index-mems
OURS=$HOME/research/rlbwt-tricks/chr19/c33; PAT=$HOME/research/rlbwt-tricks/chr19/queries.fa; RUN="taskset -c 9"
for L in 1000 100; do
  SUM=$W/summary.L$L.txt; note() { echo "$*" | tee -a "$SUM"; }
  if [ $L = 1000 ]; then cp summary.L1000.partial.txt $SUM; else
    : > $SUM; note "STPD-index ($(git -C $HOME/research/rlbwt-tricks/STPD-index log --oneline -1)) against katka ($(git -C $HOME/research/mem-finding log --oneline -1)); L = $L; patterns $PAT"
    note "their index: $(ls -l st.stpdi | awk '{printf "%.1f MiB", $5 / 1048576}') (built for L = 1000, reused)"
    for name in bml bml-locate; do opt=""; [ $name = bml-locate ] && opt=--locate
      echo "[$(date +%T)] theirs: $name"; /usr/bin/time -f "%e s wall, peak %M KB" $RUN $DRV st.stpdi $PAT $L $opt --pieces st.pieces > theirs.$name.tsv 2> theirs.$name.log
      note "theirs $name: $(grep '^stpd-index-mems' theirs.$name.log | sed 's/^stpd-index-mems: //') | $(tail -1 theirs.$name.log)"; done
  fi
  : > theirs.fb.tsv; note "theirs fb: not run (stuck for 54 min on L = 1000; see README)"
  for spec in "plain --csa" "twolevel-blocks --csa --twolevel" "stpd-nointervals --stpd --no-intervals" "stpd --stpd" "stpd-twolevel --csa --stpd --twolevel" "stpd-locate --stpd --locate" "fb-twolevel-blocks --csa --twolevel --fb"; do
    set -- $spec; name=$1; shift; echo "[$(date +%T)] ours: $name"
    /usr/bin/time -f "%e s wall, peak %M KB" $RUN $K find -L $L --format mems --stats "$@" $OURS $PAT ours.$name.tsv 2> ours.$name.log
    note "ours $name: $(grep -o 'reads=[0-9]* MEMs=[0-9]*' ours.$name.log), $(grep -o 'wall [0-9.]* s' ours.$name.log); $(grep searches ours.$name.log) | $(tail -1 ours.$name.log)"
  done
  python3 - "$SUM" < $HOME/research/rlbwt-tricks/chr19/stpd33_compare.py | tee -a "$SUM"
  mkdir -p L$L; cp theirs.*.tsv theirs.*.log ours.*.tsv ours.*.log L$L/
done; echo "== done $(date +%T)"
