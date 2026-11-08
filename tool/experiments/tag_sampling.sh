#!/bin/bash
# Backward steps on the run-length BWT, without the RLCSA, with the full tag array and with tags sampled at several s, against
# the explicit RLCSA with the full tag array: the protocol of results/silva/timing_repeats.txt, one thread, pinned, on all the
# reads, REPEATS times (tool/results/silva/README.md, "RLBWT with sampled tags").  The configurations take turns within each
# repeat, and every output is compared with the RLCSA's.
#
#   tag_sampling.sh OUT [s ...]        (default s: 1 2 4 8; 1 is the full tag array, IDX.tag)
#
# Environment: IDX (~/research/katka/work/bac), L (30), F (10), READS (katka/silva/work/compare/reads_every50), REPEATS (5), PIN (9), MEMMAX (11G).
set -uo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd); TOP=$(cd "$HERE/.." && pwd); K=$HERE/bin/katka
OUT=$1; shift; SS=${*:-1 2 4 8}; mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
DG=${IDX:-$HOME/research/katka/work/bac}; L=${L:-30}; F=${F:-10}; READS=${READS:-$TOP/katka/silva/work/compare/reads_every50}; REPEATS=${REPEATS:-5}; PIN=${PIN:-9}
REGIONS="V1_V2 V3_V4 V4_V4 V4_V5"
CAP="systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0"
for s in $SS; do
  [ "$s" = 1 ] && continue
  [ -s "$DG.s$s.tag" ] || { echo "building s = $s"; $CAP "$K" index --tags --sample "$s" --no-rmq --tag-out "$DG.s$s.tag" "$DG" > "$OUT/build.s$s.log" 2>&1 || { echo "FAILED: build s = $s"; exit 1; }; }
done
one() { # name rep find-options...
  local name=$1 rep=$2; shift 2
  local d=$OUT/$name.r$rep; mkdir -p "$d"
  for r in $REGIONS; do for m in 1 2; do
    $CAP taskset -c "$PIN" /usr/bin/time -f "peak %M KB" "$K" find -L "$L" -c -F "$F" "$@" -j 1 "$DG" "$READS/${r}_$m.fq" "$d/${r}_$m.out" > "$d/${r}_$m.log" 2>&1 \
      || { echo "FAILED: $name $r $m"; exit 1; }
    if [ "$name" = csa ]; then [ "$rep" = 1 ] && mv "$d/${r}_$m.out" "$OUT/ref.${r}_$m.out" || { cmp -s "$d/${r}_$m.out" "$OUT/ref.${r}_$m.out" || echo "DIFFERS: $name r$rep ${r}_$m"; rm -f "$d/${r}_$m.out"; }
    else cmp -s "$d/${r}_$m.out" "$OUT/ref.${r}_$m.out" && echo "$name r$rep ${r}_$m identical" >> "$OUT/identical.log" || echo "DIFFERS: $name r$rep ${r}_$m"; rm -f "$d/${r}_$m.out"; fi
  done; done
}
for rep in $(seq 1 "$REPEATS"); do
  echo "repeat $rep, $(date +%T)"
  one csa "$rep" --csa
  for s in $SS; do
    t=$DG.tag; [ "$s" = 1 ] || t=$DG.s$s.tag
    one "s$s" "$rep" -T "$t"
    echo "  s = $s done, $(date +%T)"
  done
done
rm -f "$OUT"/ref.*.out; echo "== done $(date +%T)"
