#!/bin/bash
# SILVA, preset silva (L = 30, -c --csa -F 10), one thread, 8,000 reads; 5 repeats of each variant (tool's wall, without loading)
cd ~/tool/silva/profile; K=~/rz/MEM-finding/tool/bin/katka; P=~/tool/silva/bac
run() { local name=$1; shift; for i in 1 2 3 4 5; do
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 $K find -j 1 "$@" $P reads.fq $name.out 2> $name.$i.err || { echo "FAILED $name"; exit 1; }
  grep -o 'wall [0-9.]* s, [0-9.]* us/read' $name.$i.err | sed "s/^/$name $i: /"; done; }
echo "== $(date +%T)"
run full   --preset silva
run notags -L 30 --csa -F 10 --no-tags
run noio   -L 1000 -c --csa -F 10
run mems   -L 30 --csa -F 10 --format mems
systemd-run --user --scope -q -p MemoryMax=11G $K find -j 1 --preset silva --stats $P reads.fq stats.out 2> stats.err; tail -2 stats.err
cmp -s full.out stats.out && echo "stats run: same output"
echo "== done $(date +%T)"
