#!/bin/bash
# Regression for the STPD sample and the exported index files against the paper's code
# (https://github.com/regindex/STPD-index): from the BWT runs of S = rev(T), katka must produce
#   - the sample: their st-colex- sampling of T (stpd-sampling -c), position for position, except for their
#     position 0 (which katka's own sample leaves out, and stpd-export adds);
#   - T.colex_m, T.rbwt, T.pa and T.lcs: byte for byte what stpd-sampling -c -P writes (the files stpd-build
#     reads to build its index).
#   regress_stpd.sh STPDBIN TEXT WORKDIR
# STPDBIN holds STPD-index's stpd-sampling (its build/CLI); TEXT is a file over A, C, G, T only (e.g. its
# data/yeast.txt), whose newlines are dropped.
set -euo pipefail
OLD=$(cd "$1" && pwd); TEXT=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); W=$3
HERE=$(cd "$(dirname "$0")/.." && pwd)
rm -rf "$W"; mkdir -p "$W"; cd "$W"
tr -d '\n' < "$TEXT" > T.txt
"$OLD/stpd-sampling" -i T.txt -o T.colex -c > /dev/null 2>&1
cp T.txt theirs.txt; "$OLD/stpd-sampling" -i theirs.txt -o theirs.txt.colex_m -c -P > /dev/null 2>&1
cp T.txt r.T; python3 -c "open('r.S','w').write(open('r.T').read()[::-1])"
"$HERE/bin/katka" bwt r 2> /dev/null
"$HERE/bin/katka" index --stpd r 2>&1 | grep STPD
"$HERE/bin/katka-test" stpd-vs r T.colex
"$HERE/bin/katka" stpd-export r 2> /dev/null
fails=0
for e in colex_m rbwt pa lcs; do
  if cmp -s "r.T.$e" "theirs.txt.$e"; then echo "  same: T.$e"; else echo "  DIFFERENT: T.$e"; fails=$((fails + 1)); fi
done
[ $fails = 0 ] || exit 1
