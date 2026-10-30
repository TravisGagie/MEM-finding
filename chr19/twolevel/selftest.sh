#!/bin/bash
# Builds everything and checks, on a small synthetic collection, that two-level search
# with and without block jumps, with either continuation, finds the same MEMs, with the
# same intervals, as plain backward search.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
B="$HERE/build"
cmake -S "$HERE" -B "$B" -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build "$B" -j4 --target build_twolevel query_twolevel count_strings build_index > /dev/null
T="$HERE/selftest"; rm -rf "$T"; mkdir -p "$T"
(cd "$T" && python3 "$HERE/gen_test.py" 3 200000 0.0005 > /dev/null)
"$B/rs_build/build_index" --text "$T/collection.txt" --reference-length 200000 --out "$T/idx" --ell 25 --theta 0 > /dev/null 2>&1
"$B/rs_build/build_index" --text "$T/collection.txt" --reference-length 200000 --out "$T/idxr" --ell 25 --theta 0 --reversed > /dev/null 2>&1
for c in smallest frequent; do
  "$B/build_twolevel" --text "$T/collection.txt" --rlbwt "$T/idx/rlbwt.bin" --reversed-rlbwt "$T/idxr/rlbwt.bin" \
    --out "$T/lv_$c" --continuation "$c" 2> "$T/build_$c.log"
  "$B/query_twolevel" --index "$T/idx" --reversed-index "$T/idxr" --level "$T/lv_$c.text.lvl" --reversed-level "$T/lv_$c.rev.lvl" \
    --patterns "$T/queries.txt" --shortcuts ell25_theta0 --min-length 1,100,1000 --repeat 1 > "$T/query_$c.txt" 2> "$T/query_$c.log" \
    || { echo "self-test FAILED ($c; see $T/query_$c.log)"; exit 1; }
  if grep -q MISMATCH "$T/query_$c.log"; then echo "self-test FAILED ($c)"; exit 1; fi
done
echo "self-test passed: the same MEMs and intervals as plain search"
