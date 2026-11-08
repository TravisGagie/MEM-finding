#!/bin/bash
# Another held-out split, as results/silva/heldout/ (heldout.py, then heldout2.sh) made the first: for each aquatic genus a
# random 10% of its sequences (at least 1) held out with seed SEED, KATKA's default index on the rest (timed, each step
# alone), reads simulated from the held-out sequences only, every credit scheme, and the same reads against the full index.
#   heldout_split.sh N SEED       -> ~/research/katka/work/heldoutN; the runs (10 M pairs per region at the peak) on the external
#                                    drive, linked as heldoutN/runs; the sampled reads in heldout_samplesN on the drive
set -uo pipefail
N=$1; SEED=$2; T=$(cd "$(dirname "$0")/../.." && pwd); H=$HOME/research/katka/work/heldout$N; X=$HOME/ext/moved_2026-11-05
K=$T/bin/katka; BB=$T/../katka/deps/Big-BWT/bigbwt; M=$T/../katka; S=$M/silva; Q=$HOME/research/katka/work/queue
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"; TM="/usr/bin/time -f %e_%U_%S_%M -a -o $H/build.time"
mkdir -p $H; : > $H/build.time
sed -e "s|~/research/katka/work/heldout\"|~/research/katka/work/heldout$N\"|" -e "s|random.Random(20261107)|random.Random($SEED)|" $T/results/silva/heldout/heldout.py > $H/heldout.py
grep -q "heldout$N\"" $H/heldout.py && grep -q "Random($SEED)" $H/heldout.py || { echo "could not set the split"; exit 1; }
python3 $H/heldout.py || exit 1
cd $H
echo prep >> build.time; $CAP $TM $K prep list.tsv bac > prep.log 2>&1 || exit 1
echo bigbwt >> build.time; $CAP $TM $BB bac.S > bigbwt.log 2>&1 || exit 1; mv bac.S.bwt bac.bwt; rm -f bac.S.log
echo index >> build.time; $CAP $TM $K index --preset silva bac > index.log 2>&1 || exit 1
rm -f bac.S bac.bwt bac.tagruns          # intermediates, as large as the rest together
# the genera with usable held-out amplicons, per region (as heldout2.sh)
$CAP python3 $T/scripts/sparse_genera.py --mm ~/research/datasets/MicrobeMixer --fasta $H/heldout.fasta --tax $S/exp1_data/tax_slv_ssu_138.1.txt --doc-to-traversal $S/work/ref/rna/doc_to_traversal.txt --nseq $H/nseq.txt --primers-dir ~/research/datasets/MicrobeMixer/data --out $H/genera.tsv > $H/genera.log 2>&1 || exit 1
H=$H python3 - <<'PY'
import os
H=os.environ['H']; S=os.path.expanduser('~/research/mem-finding/katka/silva')
rows={l.split('\t')[2]:l.rstrip('\n').split('\t') for l in list(open(H+'/genera.tsv'))[1:]}
ab=list(open(S+'/work/abundance/aquatic.tsv')); head,ab=ab[0],[l.rstrip('\n').split('\t') for l in ab[1:]]
R=["V1_V2","V3_V4","V4_V4","V4_V5"]; os.makedirs(H+'/ab',exist_ok=True)
with open(H+'/dropped.txt','w') as d:
    for j,r in enumerate(R):
        keep=[]
        for a in ab:
            x=rows[a[0]]; nreg,nlong=int(x[3+2*j]),int(x[4+2*j])
            if nreg>0 and nlong/nreg>0.55: keep.append(a)
            else: d.write("%s\t%s\t%d regions, %d usable\n"%(r,a[0],nreg,nlong))
        open(H+'/ab/%s.tsv'%r,'w').write(head+''.join('\t'.join(a)+'\n' for a in keep))
        print(r,len(keep),'of',len(ab),'genera kept')
PY
mkdir -p $X/heldout${N}_runs; ln -sfn $X/heldout${N}_runs $H/runs
cd $T; AB="$H/ab/%r.tsv" IDX=$H/bac NSEQ=$H/nseq.txt FASTA=$H/heldout.fasta SAMPLES=$X/heldout_samples$N experiments/aquatic_art_profiles.sh $H/runs "MSv1.1" || exit 1
Y=$X/heldout_samples$N/MSv1.1; O=$H/full; rm -rf $O; mkdir -p $O
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  for m in 1 2; do $CAP $K find --preset silva -j 12 ~/research/katka/work/bac $Y/${r}_$m.fq.gz $O/$m.out 2> /dev/null; done
  for spec in lca:lca: one:one: cliffs:cliffs: list:list: norm:norm: freq:freq: prior:prior: mix:mix: em:em: bayesflat:bayes:--prior+flat+--epsilon+1e-6,1e-4,1e-5,1e-7,1e-8 bayesauto:bayesauto: bayesem:bayesem:--em-iters+20; do
    name=${spec%%:*}; rest=${spec#*:}; mode=${rest%%:*}; opts=${rest#*:}; opts=${opts//+/ }
    $CAP python3 scripts/score_silva.py --mode $mode $opts --levels genus --nseq $Q/nseq.txt --mate1-listings $O/1.out --mate2-listings $O/2.out --doc-id-to-traversal $S/work/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt --readset-truthset $H/runs/MSv1.1/$r/${r}_seqtax.txt --trav-to-length $S/work/ref/trav_to_length.txt --exp-src $M/cliffy-experiments/src --output $O/$name.$r.csv --threads 12 > $O/$name.$r.csv.log 2>&1 || echo "FAILED $name $r"
  done; echo "full index $r $(date +%T)"
done; rm -f $O/1.out $O/2.out; echo "== heldout$N done $(date +%T)"
