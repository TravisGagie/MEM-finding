#!/bin/bash
P=~/research/katka/work/pc; D=$P/movi; cd $P; M=~/research/deps/movi_build/bin/movi-regular; K=~/research/mem-finding/tool/bin/katka; R=~/research/mem-finding/katka/silva/work/compare/reads_every50; BR=$HOME/research/deps/bigrepair/bigrepair
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
for s in 2 4 8 16; do [ -s bac.s$s.tag ] || $CAP $K index --tags --sample $s --no-rmq --tag-out bac.s$s.tag bac > s$s.log 2>&1; echo "s=$s $? $(stat -c %s bac.s$s.tag)"; done
date +"== grammar %T"; [ -s g.tag ] || $CAP /usr/bin/time -f "built-in: %e s, peak %M KB" $K index --tags --grammar --no-rmq --tag-out g.tag bac > g.log 2>&1; echo "g $? $(tail -1 g.log | cut -c1-60)"
python3 ~/research/mem-finding/tool/scripts/tagruns_to_int32.py bac.tagruns tags.int32 > /dev/null 2>&1; timeout 2h $CAP $BR -i -r -m 0 tags.int32 > bg0.log 2>&1; echo "bigrepair $?"; $CAP $K index --tags --no-rmq --bigrepair tags.int32 --tag-out bg.tag bac > bg.log 2>&1; echo "bg $? $(stat -c %s bg.tag)"; rm -f tags.int32*
date +"== timing %T"; cd $D; mkdir -p t; export MV_PLAIN=1; TMP=/tmp/mv.out
one() { name=$1; tag=$2; rep=$3; d=t/$name.r$rep; mkdir -p $d
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do for m in 1 2; do f=${r}_$m
    if [ "$name" = ours ]; then $CAP taskset -c 9 /usr/bin/time -f "peak %M KB" $K find -L 30 -c --csa -F 10 -j 1 $P/bac $R/$f.fq $TMP > $d/$f.log 2>&1
    else KATKA_OUT=$TMP KATKA_TAG=$tag $CAP taskset -c 9 /usr/bin/time -f "peak %M KB" $M query -i idx -r split/$f.fa --mem --min-mem-length 30 --ftab-k 10 -t 1 -o /tmp/mvx > $d/$f.log 2>&1; fi
    if [ $rep = 1 ]; then if [ "$name" = ours ]; then cp $TMP ref.$f; else python3 ~/research/katka/work/mcmp3.py ref.$f $TMP $name $f; fi; fi
  done; done; }
for rep in 1 2 3 4 5; do echo "repeat $rep $(date +%T)"
  one ours "" $rep
  for c in full:$P/bac.tag s2:$P/bac.s2.tag s4:$P/bac.s4.tag s8:$P/bac.s8.tag s16:$P/bac.s16.tag grammar:$P/g.tag bigrepair:$P/bg.tag; do one ${c%%:*} ${c#*:} $rep; echo "  ${c%%:*} $(date +%T)"; done
done; rm -f ref.* $TMP; echo "== done $(date +%T)"
