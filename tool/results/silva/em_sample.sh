#!/bin/bash
cd ~/research/mem-finding/tool; Q=~/research/katka/work/queue; M=$PWD/../katka; S=$M/silva; K=bin/katka; X=~/ext/moved_2026-11-05; T=~/research/katka/work/em_tmp; mkdir -p $T ~/research/katka/work/emsample
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
sc() { L=$1; m1=$2; m2=$3; tr=$4; o=~/research/katka/work/emsample/$L.csv; [ -s $o ] && return
  $CAP python3 scripts/score_silva.py --mode bayesem --em-iters 20 --epsilon 1e-6 --em-sample 1000,10000,100000 --em-smooth 0,p0.5,p1,m0.01,m0.1 --em-seeds 1,2,3 --mate1-listings $m1 --mate2-listings $m2 --doc-id-to-traversal $S/work/ref/rna/doc_to_traversal.txt \
    --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt --readset-truthset $tr --trav-to-length $S/work/ref/trav_to_length.txt --exp-src $M/cliffy-experiments/src --output $o --threads 12 --levels genus --nseq $Q/nseq.txt > $o.log 2>&1 || echo "FAILED $L"; echo "$(date +%T) $L"; }
regen() { for m in 1 2; do $CAP $K find --preset silva -j 12 ~/research/katka/work/bac "$1/${2}_$m.fq.gz" $T/$m.out 2> /dev/null || echo "FIND FAILED $1 $2"; done; }
for r in V1_V2 V3_V4 V4_V4 V4_V5; do sc aquatic.1.$r $Q/timing/default.${r}_1.out $Q/timing/default.${r}_2.out $S/work/reads/aquatic/$r/${r}_seqtax.txt; done
for set in uniform small; do for k in 1 2 3; do for r in V1_V2 V3_V4 V4_V4 V4_V5; do regen $X/sparse_samples/$set/MSv1.$k $r; sc $set.$k.$r $T/1.out $T/2.out ~/research/katka/work/sparse/$set/MSv1.$k/$r/${r}_seqtax.txt; done; done; done
rm -rf $T; echo "== em done $(date +%T)"
