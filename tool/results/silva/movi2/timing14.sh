#!/bin/bash
# Item 4: KATKA's frontier configurations, Kraken 2 and Tagger in one rotation, 5 repeats.  Every tool timed the same way:
# wall time (GNU time) on all the reads minus wall time on a single read (loading and setup), 1 thread pinned to CPU 9,
# then 12 threads.  Logs: ~/research/katka/work/t14/<config>.p<threads>.r<rep>/<file>.{all,one}.{log,time}
P=~/research/katka/work/pc; D=$P/movi; MV=~/research/deps/movi_build/bin/movi-regular; K=~/research/mem-finding/tool/bin/katka
R=~/research/mem-finding/katka/silva/work/compare/reads_every50; TG=~/research/deps/tagger/build/tagger; O=~/research/katka/work/t14; mkdir -p $O /tmp/one
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"; TM="/usr/bin/time -f %e_%U_%S_%M"
F="V1_V2_1 V1_V2_2 V3_V4_1 V3_V4_2 V4_V4_1 V4_V4_2 V4_V5_1 V4_V5_2"
for f in $F; do head -4 $R/$f.fq > /tmp/one/$f.fq; head -2 $D/split/$f.fa > /tmp/one/$f.fa; done
export MV_PLAIN=1
run() {   # config kind arg threads rep
  local n=$1 k=$2 a=$3 t=$4 rep=$5 d=$O/$1.p$4.r$5 pin=""; mkdir -p $d; [ $t = 1 ] && pin="taskset -c 9"
  local need=16; [ $k = kraken2 ] && need=8; [ "$(find $d -name '*.time' -size +0 | wc -l)" -ge $need ] && return 0   # done already
  for w in one all; do
    if [ $k = kraken2 ]; then for r in V1_V2 V3_V4 V4_V4 V4_V5; do
        i1=$R/${r}_1.fq; i2=$R/${r}_2.fq; [ $w = one ] && { i1=/tmp/one/${r}_1.fq; i2=/tmp/one/${r}_2.fq; }
        $CAP $pin $TM -o $d/$r.$w.time kraken2 --db $a --threads $t --paired $i1 $i2 > /tmp/t14.out 2> $d/$r.$w.log; done; continue; fi
    for f in $F; do
      i=$R/$f.fq; [ $w = one ] && i=/tmp/one/$f.fq
      case $k in
        csa|efcsa) $CAP $pin $TM -o $d/$f.$w.time $K find -L 30 -c --$k -F 10 -T $a -j $t $P/bac $i /tmp/t14.out > $d/$f.$w.log 2>&1 ;;
        movi) m=$D/split/$f.fa; [ $w = one ] && m=/tmp/one/$f.fa
              (cd $D && KATKA_OUT=/tmp/t14.out KATKA_TAG=$a $CAP $pin $TM -o $d/$f.$w.time $MV query -i idx -r $m --mem --min-mem-length 30 --ftab-k 10 -t $t -o /tmp/mvx > $d/$f.$w.log 2>&1) ;;
        tagger) Rr=""; tt=$t; [ $t -gt 1 ] && Rr="-R"; [ $w = one ] && { tt=1; Rr=""; }   # Tagger with several threads hangs on a one-read file
              $CAP $pin $TM -o $d/$f.$w.time $TG -f $i -r ~/research/katka/work/tagger/full/silva -L $a -t $tt $Rr -o /tmp/t14.out > $d/$f.$w.log 2>&1 ;;
      esac
      if [ $w = all ] && [ $rep = 1 ] && [ $t = 1 ] && [ $k != tagger ] && [ $k != movi ]; then
        [ -s $O/ref.$f ] || $CAP $K find --preset silva -j 12 $P/bac $R/$f.fq $O/ref.$f > /dev/null 2>&1
        cmp -s /tmp/t14.out $O/ref.$f && echo "$n $f: identical to find --preset silva" || echo "$n $f DIFFERS"; fi
    done
  done; }
C="csa_full:csa:$P/bac.tag csa_s2:csa:$P/bac.s2.tag csa_grammar:csa:$P/g.tag csa_s4:csa:$P/bac.s4.tag csa_s8:csa:$P/bac.s8.tag ef_grammar:efcsa:$P/g.tag ef_s4:efcsa:$P/bac.s4.tag ef_s8:efcsa:$P/bac.s8.tag movi_full:movi:$P/bac.tag movi_grammar:movi:$P/g.tag kraken2:kraken2:$HOME/research/katka/work/k2/std tagger_L25:tagger:25 tagger_L30:tagger:30"
for rep in 1 2 3 4 5; do echo "repeat $rep $(date +%T)"
  for t in 1 12; do for c in $C; do IFS=: read n k a <<< "$c"; run $n $k $a $t $rep; done; done; echo "  done $(date +%T)"
done; rm -f /tmp/t14.out; echo "== done $(date +%T)"
