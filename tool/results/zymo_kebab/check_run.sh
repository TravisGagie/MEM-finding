#!/bin/bash
# KeBaB revision on Zymo (tool/NOTES.md, "KeBaB: --kebab-check T t and --kebab-top t --fixed-L"), tool at dc63572+
cd ~/zymo/kebab2; K=~/research/mem-finding/tool/bin/katka; H=~/zymo/hpc/zy; R=~/zymo/all.fq; Z="-U -D 14 -c --csa -j 12 --stats"
S=~/research/mem-finding/tool/scripts; KR=~/research/mem-finding/katka-readable/scripts; T=~/zymo/all.truth; B=~/zymo/hpc/zy.tbl; P1=~/zymo/kebab
step() { local name=$1; shift; echo "== $name $(date +%T)"
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "$name: %e s wall, %U s user, peak %M KB" "$@" > $name.err 2>&1
  local rc=$?; grep -E 'reads=|kebab|wall,' $name.err; echo "$name exit $rc"; [ $rc = 0 ] || { echo "== STOPPED at $name"; exit 1; }; }
for t in 1 3 5; do TT=$((2 * t))
  step e.$TT.$t  $K find $Z -L 20 -K 20 --kebab-check $TT $t          --read-stats e.$TT.$t.tsv  $H $R e.$TT.$t
  step es.$TT.$t $K find $Z -L 20 -K 20 --kebab-check $TT $t --strict --read-stats es.$TT.$t.tsv $H $R es.$TT.$t
  step f.$t      $K find $Z -L 25 -K 20 --kebab-top $t --fixed-L $H $R f.$t
done
step g.L25 $K find $Z -L 25 -K 20 $H $R g.L25
echo "== scoring and comparisons $(date +%T)"
for t in 1 3 5; do TT=$((2 * t))
  for r in e.$TT.$t es.$TT.$t; do
    echo "-- $r: maxmem_roc (longest MEMs vote)"; python3 -I $S/maxmem_roc_stream.py $r $T $B $R | sed -n '3p'
    for ref in c d; do echo "-- $r against $ref.$t"; python3 -I $S/compare_longest_t.py $r <(zstdcat $P1/$ref.$t.zst) $R $t; done
  done
done
for r in f.1 f.3 f.5 g.L25; do
  echo "-- $r: every MEM votes"; python3 -I $KR/credit_schemes.py $r $T $B | head -3; python3 -I ~/zymo/score_spu2.py $r $T $B | head -1
done
echo "== compressing $(date +%T)"; for f in e.* es.* f.[135] g.L25; do case $f in *.err) ;; *) zstd -q -T6 -19 --rm $f;; esac; done; ls -l; df -h / | tail -1
echo "== done $(date +%T)"
