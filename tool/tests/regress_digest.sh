#!/bin/bash
# Regression for minimizer digests against KATKA's original code (katka/rz-index: rz-prep -B,
# rz-classify -B with a tag array): the same digested text, table and byte map, and the same listings,
# byte for byte, on the same BWT.
#   regress_digest.sh RZBIN DATADIR WORKDIR
# RZBIN holds rz-prep, rz-build, rz-tagbuild, rz-csabuild and rz-classify (make SDSL=... in a copy of
# katka/rz-index); DATADIR has list.tsv and reads.fq.
set -euo pipefail
OLD=$(cd "$1" && pwd); DATA=$(cd "$2" && pwd); W=$3
NEW=$(cd "$(dirname "$0")/.." && pwd)/bin/katka
rm -rf "$W"; mkdir -p "$W"; cd "$W"
fails=0
cmpout() { if cmp -s "$2" "$3"; then echo "  same: $1"; else echo "  DIFFERENT: $1"; fails=$((fails + 1)); fi; }
(cd "$DATA" && "$OLD/rz-prep" -s -B 4,11 list.tsv "$W/o") > /dev/null 2>&1
(cd "$DATA" && "$NEW" prep -s --digest 4,11 list.tsv "$W/n") 2> /dev/null
cmpout "digested S" o.S n.S; cmpout "tbl" o.tbl n.tbl; cmpout "byte map" o.map n.map
"$NEW" bwt n 2> /dev/null; cp n.bwt o.bwt
"$OLD/rz-build" -G o.S o.bwt - - o.tbl o > /dev/null 2>&1
"$OLD/rz-build" -A o.S o.bwt - - o.tbl o > /dev/null 2>&1
"$OLD/rz-tagbuild" o 1 o.tag > /dev/null 2>&1
"$OLD/rz-csabuild" o.bwt o.csa > /dev/null 2>&1
"$NEW" index --csa --tags n 2> /dev/null
for L in 11 15 20 30; do
  for F in 0 2; do
    rm -f old.out new.out
    FO=(); [ $F != 0 ] && FO=(-F $F)
    RZ_COUNTS=1 "$OLD/rz-classify" -L $L -l -T o.tag -C o.csa -B o.map "${FO[@]}" o.rz "$DATA/reads.fq" old.out 2> /dev/null
    "$NEW" find -L $L -c --csa "${FO[@]}" n "$DATA/reads.fq" new.out 2> /dev/null
    cmpout "L=$L counts, -F $F" old.out new.out
  done
  rm -f old.out new.out
  "$OLD/rz-classify" -L $L -l -T o.tag -B o.map o.rz "$DATA/reads.fq" old.out 2> /dev/null
  "$NEW" find -L $L -l n "$DATA/reads.fq" new.out 2> /dev/null
  cmpout "L=$L list, RLBWT" old.out new.out
done
if [ $fails = 0 ]; then echo "regression against rz-index digests: passed"; else echo "regression against rz-index digests: $fails differences"; exit 1; fi
