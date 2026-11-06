#!/bin/bash
# Using the lookup table at every step, against the current search, on the undigested SILVA index, by k: (a) find --absent, (b) find
# --narrow, and both, with the compact table (--compact), and at PLAIN (10) with the plain table too.  The protocol of
# results/silva/timing_repeats.txt, one thread, pinned, on all the reads, REPEATS times, the configurations taking turns within
# each repeat.  Every output is compared with the current search's at -F 10 (find --preset silva).  Then one run of the current
# search and of (a) and (b) together with the Elias-Fano RLCSA, and the probes per backward step and the backward steps, if PROBES
# names a binary built with -DKATKA_COUNT_PROBES.  The tables' sizes are in the logs.
#
#   narrow_timing.sh OUT [k ...]        (default k: 8 10 11 12)
#
# Environment: IDX (~/tool/silva/bac), READS (katka/silva/work/compare/reads_every50), REPEATS (5), PIN (9), MEMMAX (11G), PROBES.
set -uo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd); TOP=$(cd "$HERE/.." && pwd); K=$HERE/bin/katka
OUT=$1; shift; KS=${*:-8 10 11 12}; mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
IDX=${IDX:-$HOME/tool/silva/bac}; READS=${READS:-$TOP/katka/silva/work/compare/reads_every50}; REPEATS=${REPEATS:-5}; PIN=${PIN:-9}
REGIONS="V1_V2 V3_V4 V4_V4 V4_V5"
CAP="systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0"
one() { # binary name rep find-options...
  local bin=$1 name=$2 rep=$3; shift 3
  local d=$OUT/$name.r$rep; mkdir -p "$d"
  for r in $REGIONS; do for m in 1 2; do
    $CAP taskset -c "$PIN" /usr/bin/time -f "peak %M KB" "$bin" find -L 30 -c "$@" -j 1 "$IDX" "$READS/${r}_$m.fq" "$d/${r}_$m.out" > "$d/${r}_$m.log" 2>&1 \
      || { echo "FAILED: $name $r $m"; exit 1; }
    if [ ! -s "$OUT/ref.${r}_$m.out" ]; then mv "$d/${r}_$m.out" "$OUT/ref.${r}_$m.out"
    else cmp -s "$d/${r}_$m.out" "$OUT/ref.${r}_$m.out" && echo "$name r$rep ${r}_$m identical" >> "$OUT/identical.log" || echo "DIFFERS: $name r$rep ${r}_$m"; rm -f "$d/${r}_$m.out"; fi
  done; done
}
one "$K" ref 1 --csa -F 10
C="--csa --compact"
for rep in $(seq 1 "$REPEATS"); do
  echo "repeat $rep, $(date +%T)"
  for k in $KS; do
    one "$K" cur$k "$rep" --csa -F "$k"
    one "$K" c$k "$rep" $C -F "$k"
    one "$K" ca$k "$rep" $C -F "$k" --absent
    one "$K" cb$k "$rep" $C -F "$k" --narrow
    one "$K" cab$k "$rep" $C -F "$k" --absent --narrow
    if [ "$k" = "${PLAIN:-10}" ]; then
      one "$K" pa$k "$rep" --csa -F "$k" --absent
      one "$K" pb$k "$rep" --csa -F "$k" --narrow
      one "$K" pab$k "$rep" --csa -F "$k" --absent --narrow
    fi
    echo "  k = $k done, $(date +%T)"
  done
done
for k in $KS; do
  one "$K" efcur$k 1 --efcsa -F "$k"
  one "$K" efcab$k 1 --efcsa --compact -F "$k" --absent --narrow
  echo "  Elias-Fano, k = $k done, $(date +%T)"
done
if [ -n "${PROBES:-}" ]; then
  for k in $KS; do
    one "$PROBES" pcur$k 1 --csa -F "$k"
    one "$PROBES" pca$k 1 $C -F "$k" --absent
    one "$PROBES" pcb$k 1 $C -F "$k" --narrow
    one "$PROBES" pcab$k 1 $C -F "$k" --absent --narrow
  done
  echo "probes done, $(date +%T)"
fi
rm -f "$OUT"/ref.*.out; echo "== done $(date +%T)"
