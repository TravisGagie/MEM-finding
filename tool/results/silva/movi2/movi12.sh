#!/bin/bash
P=~/tool/silva/pc; D=$P/movi; cd $D; M=~/tool/movi_build/bin/movi-regular; K=~/rz/MEM-finding/tool/bin/katka; R=~/rz/MEM-finding/katka/silva/work/compare/reads_every50; export MV_PLAIN=1; TMP=/tmp/mv.out
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
one() { name=$1; kind=$2; tag=$3; rep=$4; d=t5/$name.r$rep; mkdir -p $d
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do for m in 1 2; do f=${r}_$m
    case $kind in
     csa|efcsa) $CAP taskset -c 9 /usr/bin/time -f "peak %M KB" $K find -L 30 -c --$kind -F 10 -T $tag -j 1 $P/bac $R/$f.fq $TMP > $d/$f.log 2>&1 ;;
     movi) KATKA_OUT=$TMP KATKA_TAG=$tag $CAP taskset -c 9 /usr/bin/time -f "peak %M KB" $M query -i idx -r split/$f.fa --mem --min-mem-length 30 --ftab-k 10 -t 1 -o /tmp/mvx > $d/$f.log 2>&1 ;;
    esac
    if [ $rep = 1 ]; then
      if [ $name = csa_full ]; then $CAP $K find --preset silva -j 12 $P/bac $R/$f.fq /tmp/def.out > /dev/null 2>&1; cmp -s $TMP /tmp/def.out && echo "csa_full $f: identical to find --preset silva" || echo "csa_full $f DIFFERS from the default"; cp $TMP ref.$f
      elif [ $kind = movi ]; then python3 ~/tool/silva/mcmp3.py ref.$f $TMP $name $f
      else cmp -s $TMP ref.$f && echo "$name $f: different 0" || echo "$name $f DIFFERS"; fi
    fi
  done; done; }
for rep in 1 2 3 4 5; do echo "repeat $rep $(date +%T)"
  for c in csa_full:csa:$P/bac.tag csa_grammar:csa:$P/g.tag csa_s4:csa:$P/bac.s4.tag ef_full:efcsa:$P/bac.tag ef_grammar:efcsa:$P/g.tag ef_s4:efcsa:$P/bac.s4.tag movi_full:movi:$P/bac.tag movi_grammar:movi:$P/g.tag movi_s4:movi:$P/bac.s4.tag; do
    IFS=: read n k t <<< "$c"; one $n $k $t $rep; done; echo "  done $(date +%T)"
done; rm -f ref.* $TMP /tmp/def.out; echo "== done $(date +%T)"
