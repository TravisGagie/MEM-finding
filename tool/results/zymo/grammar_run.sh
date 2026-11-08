#!/bin/bash
# Grammar-compressed per-species tags on Zymo, with bigrepair -m 0 (tool/results/zymo/README.md).  The built-in RePair isn't run:
# it needs about 25 bytes per run, 11.4 GB here, over the 11 GB cap.  From experiments/grammar_tags.sh's steps 2 to 4.
H=$HOME/zymo/hpc; P=$H/zy; W=$H/grammar; mkdir -p "$W"; cd ~/research/mem-finding
K=$PWD/tool/bin/katka; BR=$HOME/research/deps/bigrepair/bigrepair; BWT="$HOME/ext/zymo_bwt/zy.bwt"
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
step() { local n=$1; shift; echo "== $n $(date +%T)"; $CAP /usr/bin/time -f "$n: %e s wall, %U s user, peak %M KB" "$@" > "$W/$n.log" 2>&1; local rc=$?; grep -E '^(tag array|RLBWT|'"$n"':)' "$W/$n.log" | cut -c1-120; echo "$n exit $rc"; return $rc; }
echo "rlbwt before: $(sha256sum "$P.rlbwt" | cut -c1-16)"
step int32 python3 tool/scripts/tagruns_to_int32.py "$P.sp.tagruns" "$W/tags.int32" || exit 1
ls -l "$W/tags.int32" | awk '{print "tags.int32", $5}'
step bigrepair timeout 8h "$BR" -i -r -m 0 "$W/tags.int32" || exit 1
ls -l "$W/tags.int32.R" "$W/tags.int32.C" | awk '{print $9, $5}'
step index "$K" index --tags --by-tag --no-rmq --bigrepair "$W/tags.int32" --tag-out "$W/bg.tag" --bwt "$BWT" "$P" || exit 1
echo "rlbwt after: $(sha256sum "$P.rlbwt" | cut -c1-16)"
ls -l "$P.tag" "$W/bg.tag" | awk '{print $9, $5}'
head -n 400000 ~/zymo/all.fq > "$W/sample.fq"; echo "reads: $(( $(wc -l < "$W/sample.fq") / 4 ))"
step find.plain "$K" find -U -X -L 20 -D 14 -c --csa -j 12 "$P" "$W/sample.fq" "$W/out.plain"
step find.bg "$K" find -U -X -L 20 -D 14 -c --csa -T "$W/bg.tag" -j 12 "$P" "$W/sample.fq" "$W/out.bg"
cmp "$W/out.plain" "$W/out.bg" && echo "outputs byte-identical: $(wc -c < "$W/out.plain") bytes"
grep -h 'timing:' "$W/find.plain.log" "$W/find.bg.log" | cut -c1-110
rm -f "$W/out.plain" "$W/out.bg" "$W/sample.fq" "$W/tags.int32"
