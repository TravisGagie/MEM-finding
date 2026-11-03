#!/bin/bash
H=~/zymo/hpc; O=~/zymo/newrun/reg; cd $O; NEW=~/rz/MEM-finding/tool/bin/katka
step() { local name=$1; shift; echo "== $name $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$name: %e s wall, %U s user, peak %M KB" "$@" > $name.err 2>&1
  local rc=$?; tail -3 $name.err; echo "$name exit $rc"; [ $rc = 0 ] || { echo "== STOPPED at $name"; exit 1; }; }
step new.X.nokebab $NEW find -U -X -L 20 -D 14 -c --csa -j 12 $H/zy ~/zymo/all.fq new.X.nokebab
if cmp -s new.X.nokebab <(zstdcat new.X.zst); then echo "SAME as the KeBaB run: the filter loses no longest MEM"; else echo "DIFFERENT from the KeBaB run"; fi
step new.p720.L30 $NEW find -U -L 30 -D 14 -c -j 12 $H/zy ~/zymo/pfx/p720.fq new.p720.L30
python3 -I - <<'PY'
truth = dict(l.split()[:2] for l in open('/home/travis/zymo/pfx/truth'))
has = {}; name = None
for line in open('new.p720.L30'):
    if line[0] == '>': name = line[1:].strip(); has[name] = False
    elif '[' in line: has[name] = True
pos = [k for k in has if truth[k] != 'yeast']; neg = [k for k in has if truth[k] == 'yeast']
print('prefix test, 720 bases, L = 30, -D 14: %d bacterial, %d yeast prefixes' % (len(pos), len(neg)))
print('  sensitivity %.2f%%  specificity %.2f%%' % (100 * sum(has[k] for k in pos) / len(pos), 100 * sum(not has[k] for k in neg) / len(neg)))
PY
for f in new.X.nokebab new.p720.L30; do zstd -q -T6 -19 --rm $f; done; ls -l; echo "== done $(date +%T)"
