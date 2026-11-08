#!/bin/bash
# The explicit RLCSA with tags sampled at s, recovering unsampled tags by LF from the RLCSA, with no RLBWT loaded, beside the full and
# grammar-compressed tag arrays: the protocol of results/silva/timing_repeats.txt, one thread, pinned, on all the reads,
# REPEATS times, taking turns within each repeat (tool/results/silva/README.md, "RLCSA with sampled tags").  Every output is
# compared with the full tag array's.
#
#   csa_sampled_tags.sh OUT
#
# Environment: IDX (~/research/katka/work/bac; IDX.s<s>.tag for s = 2, 4, 8: see tag_sampling.sh), GRAMMAR (~/research/katka/work/queue/grammar,
# with g.tag and bg.tag), READS, REPEATS (5), PIN (9), MEMMAX (11G).
set -uo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd); TOP=$(cd "$HERE/.." && pwd); K=$HERE/bin/katka
OUT=$1; mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
IDX=${IDX:-$HOME/research/katka/work/bac}; G=${GRAMMAR:-$HOME/research/katka/work/queue/grammar}; READS=${READS:-$TOP/katka/silva/work/compare/reads_every50}; REPEATS=${REPEATS:-5}; PIN=${PIN:-9}
REGIONS="V1_V2 V3_V4 V4_V4 V4_V5"
CAP="systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0"
for rep in $(seq 1 "$REPEATS"); do
  echo "repeat $rep, $(date +%T)"
  for cfg in full:$IDX.tag s2:$IDX.s2.tag s4:$IDX.s4.tag s8:$IDX.s8.tag grammar:$G/g.tag bigrepair:$G/bg.tag; do
    name=${cfg%%:*}; tag=${cfg#*:}; d=$OUT/$name.r$rep; mkdir -p "$d"
    for r in $REGIONS; do for m in 1 2; do
      $CAP taskset -c "$PIN" /usr/bin/time -f "peak %M KB" "$K" find --preset silva -T "$tag" -j 1 "$IDX" "$READS/${r}_$m.fq" "$d/${r}_$m.out" > "$d/${r}_$m.log" 2>&1 \
        || { echo "FAILED: $name $r $m"; exit 1; }
      if [ ! -s "$OUT/ref.${r}_$m.out" ]; then mv "$d/${r}_$m.out" "$OUT/ref.${r}_$m.out"
      else cmp -s "$d/${r}_$m.out" "$OUT/ref.${r}_$m.out" && echo "$name r$rep ${r}_$m identical" >> "$OUT/identical.log" || echo "DIFFERS: $name r$rep ${r}_$m"; rm -f "$d/${r}_$m.out"; fi
    done; done
    echo "  $name done, $(date +%T)"
  done
done
rm -f "$OUT"/ref.*.out; echo "== done $(date +%T)"
