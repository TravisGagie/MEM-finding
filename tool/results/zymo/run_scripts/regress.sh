#!/bin/bash
# Full-scale regression on Zymo: KATKA (katka-readable's katka-classify on KATKA's laptop index files)
# against the new tool (katka find on the new index), same reads, same options, 12 threads, one at a time.
H=~/zymo/hpc; R=~/zymo/all.fq; O=~/zymo/newrun/reg; cd $O
OLD=~/rz/MEM-finding/katka-readable/bin/katka-classify; NEW=~/rz/MEM-finding/tool/bin/katka
step() { local name=$1; shift; echo "== $name $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$name: %e s wall, %U s user, peak %M KB" "$@" > $name.err 2>&1
  local rc=$?; tail -4 $name.err; echo "$name exit $rc"; [ $rc = 0 ] || { echo "== STOPPED at $name"; exit 1; }; }
same() { if cmp -s $1 $2; then echo "SAME: $1 $2 ($(stat -c %s $1) bytes)"; else echo "DIFFERENT: $1 $2"; cmp $1 $2 | head -2; fi; }
step old.L30 $OLD -U -L 30 -D 14 -c -C $H/katka/zy.csa -T $H/katka/zy.s1.tag -j 12 $H/katka/zy.rz $R old.L30
step new.L30 $NEW find -U -L 30 -D 14 -c --csa -j 12 $H/zy $R new.L30
same old.L30 new.L30
step new.X $NEW find --preset zymo -j 12 $H/zy $R new.X
step old.X $OLD -U -X -L 20 -K $H/katka/zy.k20.kbf -D 14 -c -C $H/katka/zy.csa -T $H/katka/zy.s1.tag -j 12 $H/katka/zy.rz $R old.X
same old.X new.X
df -h / | tail -1
echo "== compressing $(date +%T)"; for f in old.L30 new.L30 old.X new.X; do zstd -q -T6 -19 --rm $f; done; ls -l
echo "== done $(date +%T)"
