#!/bin/bash
# Regression against KATKA (katka-readable/): the same listings, byte for byte, from the old programs and
# from katka find, on the same text and BWT.
#   regress_katka.sh OLDBIN DATADIR WORKDIR
# OLDBIN holds katka-readable's katka-prep, katka-build, katka-tagbuild, katka-csabuild, katka-kbbuild and
# katka-classify (make SDSL=... in katka-readable/); DATADIR has list.tsv and reads.fq (gen_data.py, or
# KATKA's test-small with its list and reads).
set -euo pipefail
OLD=$(cd "$1" && pwd); DATA=$(cd "$2" && pwd); W=$3
NEW=$(cd "$(dirname "$0")/.." && pwd)/bin/katka
mkdir -p "$W"; cd "$W"
READS=$DATA/reads.fq
fails=0
cmpout() {   # name old new
  if cmp -s "$2" "$3"; then echo "  same: $1"; else echo "  DIFFERENT: $1"; fails=$((fails + 1)); fi
}
for hpc in 0 1; do
  U=""; [ $hpc = 1 ] && U="-U"
  p=o$hpc; q=n$hpc
  echo "== text ${U:-without homopolymer compression}"
  (cd "$DATA" && "$OLD/katka-prep" -s $U list.tsv "$W/$p") 2> /dev/null
  (cd "$DATA" && "$NEW" prep -s $U list.tsv "$W/$q") 2> /dev/null
  cmpout "S" $p.S $q.S; cmpout "tbl" $p.tbl $q.tbl
  "$NEW" bwt $q 2> /dev/null; cp $q.bwt $p.bwt          # the same BWT for both
  "$OLD/katka-build" -G -A $p.S $p.bwt $p.tbl $p > /dev/null 2>&1
  "$OLD/katka-tagbuild" $p 1 $p.tag > /dev/null 2>&1
  "$OLD/katka-tagbuild" $p 4 $p.s4.tag > /dev/null 2>&1
  "$OLD/katka-tagbuild" $p 1 $p.sp.tag -S > /dev/null 2>&1
  "$OLD/katka-csabuild" $p.bwt $p.csa > /dev/null 2>&1
  "$OLD/katka-csabuild" $p.bwt $p.efcsa ef > /dev/null 2>&1
  "$OLD/katka-kbbuild" -k 20 $p.S $p.k20.kbf > /dev/null 2>&1
  "$NEW" index --csa --efcsa --tags --kebab 20 $q 2> /dev/null
  "$NEW" index --tags --sample 4 --tag-out $q.s4.tag $q 2> /dev/null
  "$NEW" index --tags --by-tag --tag-out $q.sp.tag $q 2> /dev/null
  "$NEW" index --tags --grammar --tag-out $q.g.tag $q 2> /dev/null
  cmpout "KeBaB filter" $p.k20.kbf $q.k20.kbf
  run() {   # name, old options, new options
    local name=$1; shift
    local o=() n=()
    while [ "$1" != "--" ]; do o+=("$1"); shift; done; shift
    n=("$@")
    "$OLD/katka-classify" $U "${o[@]}" $p.rlbwt "$READS" old.out 2> /dev/null
    "$NEW" find $U "${n[@]}" $q "$READS" new.out 2> /dev/null
    cmpout "$name" old.out new.out
  }
  run "L=30 count, RLCSA, -F 10"     -L 30 -c -C $p.csa -F 10 -T $p.tag      -- -L 30 -c --csa -F 10
  run "L=20 list, RLBWT"              -L 20 -l -T $p.tag                      -- -L 20 -l
  run "L=20 range, Elias-Fano RLCSA"  -L 20 -C $p.efcsa -T $p.tag             -- -L 20 --efcsa
  run "L=25 count, s=4"               -L 25 -c -C $p.csa -T $p.s4.tag         -- -L 25 -c --csa -T $q.s4.tag
  run "L=25 list, s=4"                -L 25 -l -C $p.csa -T $p.s4.tag         -- -L 25 -l --csa -T $q.s4.tag
  run "L=25 count, by tag"            -L 25 -c -C $p.csa -T $p.sp.tag         -- -L 25 -c --csa -T $q.sp.tag
  run "L=25 count, grammar"           -L 25 -c -C $p.csa -T $p.tag            -- -L 25 -c --csa -T $q.g.tag
  run "L=30 -D 14 count"              -L 30 -D 14 -c -C $p.csa -T $p.tag      -- -L 30 -D 14 -c --csa
  run "longest -X -L 1 -D 14"         -X -L 1 -D 14 -c -C $p.csa -T $p.tag    -- -X -L 1 -D 14 -c --csa
  run "longest -X -L 20 -K"           -X -L 20 -K $p.k20.kbf -c -C $p.csa -T $p.tag -- -X -L 20 -K 20 -c --csa
  run "L=20 -K (pseudo-MEMs)"         -L 20 -K $p.k20.kbf -c -C $p.csa -T $p.tag -- -L 20 -K 20 -c --csa
  run "4 threads"                     -L 20 -c -C $p.csa -T $p.tag -j 4       -- -L 20 -c --csa -j 4
done
if [ $fails = 0 ]; then echo "regression against KATKA: passed"; else echo "regression against KATKA: $fails differences"; exit 1; fi
