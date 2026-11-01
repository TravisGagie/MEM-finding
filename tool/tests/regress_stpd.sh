#!/bin/bash
# Regression for the STPD sample against the paper's code (https://github.com/regindex/STPD-index):
# our pruned sample, built from the BWT runs of S = rev(T), must be their st-colex- sampling of T
# (stpd-sampling -c), position for position, except for their position 0.
#   regress_stpd.sh STPDBIN TEXT WORKDIR
# STPDBIN holds STPD-index's stpd-sampling (its build/CLI); TEXT is a file over A, C, G, T only (e.g. its
# data/yeast.txt), whose newlines are dropped.
set -euo pipefail
OLD=$(cd "$1" && pwd); TEXT=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); W=$3
HERE=$(cd "$(dirname "$0")/.." && pwd)
rm -rf "$W"; mkdir -p "$W"; cd "$W"
tr -d '\n' < "$TEXT" > T.txt
python3 -c "open('r.S','w').write(open('T.txt').read()[::-1])"
"$OLD/stpd-sampling" -i T.txt -o T.colex -c > /dev/null 2>&1
"$HERE/bin/katka" bwt r 2> /dev/null
"$HERE/bin/katka" index --stpd r 2>&1 | grep STPD
"$HERE/bin/katka-test" stpd-vs r T.colex
