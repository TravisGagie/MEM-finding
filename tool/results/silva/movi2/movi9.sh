#!/bin/bash
P=~/tool/silva/pc; D=$P/movi; cd $D; M=~/tool/movi_build/bin/movi-regular; K=~/rz/MEM-finding/tool/bin/katka; R=~/rz/MEM-finding/katka/silva/work/compare/reads_every50; mkdir -p t3; export MV_PLAIN=1; TMP=/tmp/mv.out
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
one() { name=$1; kind=$2; tag=$3; rep=$4; d=t3/$name.r$rep; mkdir -p $d
  for r in V1_V2 V3_V4 V4_V4 V4_V5; do for m in 1 2; do f=${r}_$m
    case $kind in
     csa) $CAP taskset -c 9 /usr/bin/time -f "peak %M KB" $K find -L 30 -c --csa -F 10 -T $tag -j 1 $P/bac $R/$f.fq $TMP > $d/$f.log 2>&1 ;;
     old) KATKA_OUT=$TMP KATKA_TAG=$tag $CAP taskset -c 9 /usr/bin/time -f "peak %M KB" $M query -i idx -r split/$f.fa --mem --min-mem-length 30 --ftab-k 10 -t 1 -o /tmp/mvx > $d/$f.log 2>&1 ;;
     sel) KATKA_OUT=$TMP KATKA_TAG=$tag KATKA_SELECT=$P/bridge.bin $CAP taskset -c 9 /usr/bin/time -f "peak %M KB" $M query -i idx -r split/$f.fa --mem --min-mem-length 30 --ftab-k 10 -t 1 -o /tmp/mvx > $d/$f.log 2>&1 ;;
    esac
    if [ $rep = 1 ]; then
      if [ "$name" = csafull ]; then cp $TMP ref.$f
      elif [ $kind = csa ]; then cmp -s $TMP ref.$f && echo "$name $f: different 0" || echo "$name $f DIFFERS"
      else python3 ~/tool/silva/mcmp3.py ref.$f $TMP $name $f; fi
    fi
  done; done; }
for rep in 1 2 3 4 5; do echo "repeat $rep $(date +%T)"
  one csafull csa $P/bac.tag $rep; one oldfull old $P/bac.tag $rep; one selfull sel $P/bac.tag $rep; echo "  full $(date +%T)"
  for s in 2 4 8; do one csas$s csa $P/bac.s$s.tag $rep; one sels$s sel $P/bac.s$s.tag $rep; done; echo "  sampled $(date +%T)"
  one csag csa $P/g.tag $rep; one selg sel $P/g.tag $rep; echo "  grammar $(date +%T)"
done; rm -f ref.* $TMP; echo "== done $(date +%T)"
