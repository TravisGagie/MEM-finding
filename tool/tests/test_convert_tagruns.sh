#!/bin/bash
# KATKA's per-genome and per-species tag runs (katka-readable's katka-tagbuild -S), converted by
# scripts/convert_katka_tagruns.py, must equal katka index's own sp.tagruns byte for byte; check_tagruns
# must find no mismatch in them and must catch one flipped tag. Runs on the self-test data.
#   test_convert_tagruns.sh OLDBIN [WORKDIR]     (OLDBIN: katka-readable/bin, built with make SDSL=...)
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd); OLD=$(cd "$1" && pwd); W=${2:-$HERE/selftest/conv}
NEW=$HERE/bin/katka; DATA=$HERE/selftest/data; CHK=$HERE/bin/check_tagruns
[ -s "$DATA/list.tsv" ] || { echo "run ./selftest.sh first (it writes selftest/data)"; exit 1; }
g++ -std=c++17 -O3 -march=native -I"$HERE/src" "$HERE/tests/check_tagruns.cpp" -o "$CHK" -lz -pthread
mkdir -p "$W"; cd "$W"; fails=0
for mode in plain hpc datasets; do
  U=""; D=""; [ $mode != plain ] && U="-U"; [ $mode = datasets ] && D="-d"; p=o$mode; q=n$mode
  (cd "$DATA" && "$OLD/katka-prep" -s $D $U list.tsv "$W/$p") 2>/dev/null
  (cd "$DATA" && "$NEW" prep -s $D $U list.tsv "$W/$q") 2>/dev/null
  "$NEW" bwt $D $q 2>/dev/null; cp $q.bwt $p.bwt
  "$OLD/katka-build" -G -A $p.S $p.bwt $p.tbl $p > build.log 2>&1
  "$OLD/katka-tagbuild" $p 1 $p.sp.tag -S > /dev/null 2>&1
  "$NEW" index --tags --by-tag --tag-out $q.sp.tag $q 2>/dev/null
  python3 -I "$HERE/scripts/convert_katka_tagruns.py" $p.sp.tagruns $p.conv > /dev/null
  if cmp -s $p.conv $q.sp.tagruns; then echo "  same: $mode (converted KATKA tag runs)"; else echo "  DIFFERENT: $mode"; fails=$((fails + 1)); fi
  "$CHK" $q $p.conv 1000000 1 > chk.out || { echo "  check FAILED: $mode"; cat chk.out; fails=$((fails + 1)); }
  python3 -I -c "
import numpy as np; d = bytearray(open('$p.conv', 'rb').read()); r = int(np.frombuffer(bytes(d[16:24]), '<u8')[0])
o = 32 + 8 * r + 4 * (r // 2); d[o] ^= 1; open('$p.bad', 'wb').write(d)"
  if "$CHK" $q $p.bad 1000000 1 > chk.out; then echo "  check MISSED a flipped tag: $mode"; fails=$((fails + 1)); else echo "  caught: one flipped tag ($mode)"; fi
done
[ $fails = 0 ] && echo "tag-run conversion test passed" || { echo "$fails failures"; exit 1; }
