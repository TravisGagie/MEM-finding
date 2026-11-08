#!/bin/bash
P=~/research/katka/work/pc; D=$P/movi; cd $D; M=~/research/deps/movi_build/bin/movi-regular; K=~/research/mem-finding/tool/bin/katka; R=~/research/mem-finding/katka/silva/work/compare/reads_every50; export MV_PLAIN=1; mkdir -p t2
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
for rep in 1 2 3 4 5; do
 for c in search rows ours; do d=t2/$c.r$rep; mkdir -p $d
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do for m in 1 2; do f=${r}_$m
   case $c in
    search) $CAP taskset -c 9 $M query -i idx -r split/$f.fa --mem --min-mem-length 30 --ftab-k 10 -t 1 -o /tmp/mvy > $d/$f.log 2>&1 ;;
    rows) KATKA_OUT=/tmp/mv.out KATKA_MEMS=1 $CAP taskset -c 9 $M query -i idx -r split/$f.fa --mem --min-mem-length 30 --ftab-k 10 -t 1 -o /tmp/mvy > $d/$f.log 2>&1 ;;
    ours) $CAP taskset -c 9 $K find -L 30 --format mems --csa -F 10 -j 1 $P/bac $R/$f.fq /tmp/mv.out > $d/$f.log 2>&1 ;;
   esac
  done; done
 done
done; rm -f /tmp/mv.out /tmp/mvy*; echo "== done $(date +%T)"
