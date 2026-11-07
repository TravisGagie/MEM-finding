#!/bin/bash
while systemctl --user is-active -q movi9; do sleep 30; done
P=~/tool/silva/pc; D=$P/movi; cd $D; M=~/tool/movi_build/bin/movi-regular; K=~/rz/MEM-finding/tool/bin/katka; R=~/rz/MEM-finding/katka/silva/work/compare/reads_every50; mkdir -p t4; export MV_PLAIN=1; TMP=/tmp/mv.out
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
one() { name=$1; kind=$2; tag=$3; rep=$4; d=t3/$name.r$rep; mkdir -p $d
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do for m in 1 2; do f=${r}_$m
    KATKA_OUT=$TMP KATKA_TAG=$tag $CAP taskset -c 9 /usr/bin/time -f "peak %M KB" $M query -i idx -r split/$f.fa --mem --min-mem-length 30 --ftab-k 10 -t 1 -o /tmp/mvx > $d/$f.log 2>&1
  done; done; }
for rep in 1 2 3 4 5; do echo "repeat $rep $(date +%T)"
  for c in olds2:$P/bac.s2.tag olds4:$P/bac.s4.tag olds8:$P/bac.s8.tag olds16:$P/bac.s16.tag oldg:$P/g.tag oldbg:$P/bg.tag; do one ${c%%:*} old ${c#*:} $rep; done
done; rm -f $TMP; echo "== done $(date +%T)"
