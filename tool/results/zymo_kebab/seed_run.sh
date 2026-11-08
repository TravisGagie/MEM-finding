#!/bin/bash
# --check-seed reruns (tool/NOTES.md): t = 3 and 5, T = 2t, with and without --strict; each cmp'd with its unseeded run
cd ~/zymo/kebab2; K=~/research/mem-finding/tool/bin/katka; H=~/zymo/hpc/zy; R=~/zymo/all.fq; Z="-U -D 14 -c --csa -j 12 --stats"
for t in 3 5; do TT=$((2 * t))
  for v in "ec e" "ecs es"; do set -- $v; new=$1.$TT.$t; old=$2.$TT.$t; extra=""; [ $1 = ecs ] && extra=--strict
    echo "== $new $(date +%T)"
    systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$new: %e s wall, %U s user, peak %M KB" $K find $Z -L 20 -K 20 --kebab-check $TT $t $extra --check-seed $H $R $new > $new.err 2>&1; echo "exit $?"
    grep -E 'kebab check|kebab \(|wall,' $new.err; grep -E 'kebab check|wall,' $old.err | sed 's/^/   unseeded: /'
    if cmp -s $new <(zstdcat $old.zst); then echo "   $new: byte-identical to $old"; else echo "   $new: DIFFERS from $old"; fi
    zstd -q -T6 --rm $new
  done
done; echo "== done $(date +%T)"
