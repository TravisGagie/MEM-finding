#!/bin/bash
# every credit scheme on the uniform-abundance reads (the listings of silva_queue.sh uniform), scored as in its credit step
M=$HOME/research/mem-finding/katka; S=$M/silva; EXP=$M/cliffy-experiments/src; Q=~/research/katka/work/queue; U=$Q/uniform/work; O=$Q/uniform/credit; mkdir -p $O
HERE=~/research/archive/mem_build-worktree/tool
for r in V1_V2 V3_V4 V4_V4 V4_V5; do for m in lca one cliffs list norm freq; do
  extra=""; [ $m = norm ] && extra="--nseq $Q/nseq.txt"
  [ -s $O/$m.$r.csv ] || python3 $HERE/scripts/score_silva.py --mode $m --mate1-listings $U/$r/out_1 --mate2-listings $U/$r/out_2 --doc-id-to-traversal $S/work/ref/rna/doc_to_traversal.txt \
    --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt --readset-truthset $U/$r/${r}_seqtax.txt --trav-to-length $S/work/ref/trav_to_length.txt --exp-src $EXP --output $O/$m.$r.csv $extra --threads 6 > $O/$m.$r.csv.log 2>&1 || echo "FAILED $m $r"
done; echo "$r done $(date +%T)"; done
