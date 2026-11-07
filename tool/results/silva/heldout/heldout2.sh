#!/bin/bash
H=~/tool/silva/heldout; cd ~/rz/MEM-finding/tool; K=$PWD/bin/katka; M=$PWD/../katka; S=$M/silva; Q=~/tool/silva/queue
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
$CAP python3 scripts/sparse_genera.py --mm ~/rz/MicrobeMixer --fasta $H/heldout.fasta --tax $S/exp1_data/tax_slv_ssu_138.1.txt --doc-to-traversal $S/work/ref/rna/doc_to_traversal.txt --nseq $H/nseq.txt --primers-dir ~/rz/MicrobeMixer/data --out $H/genera.tsv > $H/genera.log 2>&1; echo "genera $?"
python3 - <<'PY'
import os
H=os.path.expanduser('~/tool/silva/heldout'); S=os.path.expanduser('~/rz/MEM-finding/katka/silva')
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
rm -rf $H/runs ~/ext/moved_2026-11-05/heldout_samples
AB="$H/ab/%r.tsv" IDX=$H/bac NSEQ=$H/nseq.txt FASTA=$H/heldout.fasta SAMPLES=~/ext/moved_2026-11-05/heldout_samples experiments/aquatic_art_profiles.sh $H/runs "MSv1.1"; echo "runs $?"
X=~/ext/moved_2026-11-05/heldout_samples/MSv1.1; O=$H/full; rm -rf $O; mkdir -p $O
for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  for m in 1 2; do $CAP $K find --preset silva -j 12 ~/tool/silva/bac $X/${r}_$m.fq.gz $O/$m.out 2> /dev/null; done
  for spec in lca:lca: one:one: cliffs:cliffs: list:list: norm:norm: freq:freq: prior:prior: mix:mix: em:em: bayesflat:bayes:--prior+flat+--epsilon+1e-6,1e-4,1e-5,1e-7,1e-8 bayesauto:bayesauto: bayesem:bayesem:--em-iters+20; do
    name=${spec%%:*}; rest=${spec#*:}; mode=${rest%%:*}; opts=${rest#*:}; opts=${opts//+/ }
    $CAP python3 scripts/score_silva.py --mode $mode $opts --levels genus --nseq $Q/nseq.txt --mate1-listings $O/1.out --mate2-listings $O/2.out --doc-id-to-traversal $S/work/ref/rna/doc_to_traversal.txt --silva-tax-ranks $S/exp1_data/tax_slv_ssu_138.1.txt --readset-truthset $H/runs/MSv1.1/$r/${r}_seqtax.txt --trav-to-length $S/work/ref/trav_to_length.txt --exp-src $M/cliffy-experiments/src --output $O/$name.$r.csv --threads 12 > $O/$name.$r.csv.log 2>&1 || echo "FAILED $name $r"
  done; echo "full index $r $(date +%T)"
done; rm -f $O/1.out $O/2.out; echo "== heldout done $(date +%T)"
