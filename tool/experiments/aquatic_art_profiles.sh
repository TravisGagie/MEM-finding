#!/bin/bash
# The aquatic reads again under a given ART error profile: MicrobeMixer as published (6409f5d) with only `-ss` changed,
# the same abundance table, regions and read counts, every 50th pair, then every credit scheme (tool/results/silva/README.md,
# "ART's error profile").  Neither ART nor MicrobeMixer is seeded, so each repeat differs.
#
#   aquatic_art_profiles.sh OUT "PROFILE.REPEAT ..."        e.g. "MSv3.1 MSv1.2 MSv3.2"
#
# OUT/<PROFILE.REPEAT>/<region>/: credit/*.csv, the 10 M pairs' seqtax.txt and stdout.log.  The sampled reads are kept
# gzipped in SAMPLES/<PROFILE.REPEAT>/ (default OUT/samples), since they can't be regenerated; full reads and listings are deleted.
# Needs about 35 GB of disk per region at the peak.  Environment: AB (the abundance table; %r stands for the region), FASTA (the reference simulated from), S (katka/silva), IDX, NSEQ, EXP, MM, REGIONS, THREADS, MEMMAX.
set -uo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd); TOP=$(cd "$HERE/.." && pwd); K=$HERE/bin/katka
OUT=$1; PR=$2; mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd); SAMPLES=${SAMPLES:-$OUT/samples}
S=${S:-$TOP/katka/silva}; IDX=${IDX:-$HOME/tool/silva/bac}; NSEQ=${NSEQ:-$HOME/tool/silva/queue/nseq.txt}
EXP=${EXP:-$TOP/katka/cliffy-experiments/src}; MM=${MM:-$HOME/rz/MicrobeMixer}; THREADS=${THREADS:-12}
REGIONS=${REGIONS:-"V1_V2 V3_V4 V4_V4 V4_V5"}
FASTA=${FASTA:-$S/exp1_data/SILVA_138.1_SSURef_NR99_tax_silva.fasta}; TAX=$S/exp1_data/tax_slv_ssu_138.1.txt
D2T=$S/work/ref/rna/doc_to_traversal.txt; T2L=$S/work/ref/trav_to_length.txt; AB=${AB:-$S/work/abundance/aquatic.tsv}
CAP="systemd-run --user --scope -q -p MemoryMax=${MEMMAX:-11G} -p MemorySwapMax=0"
stamp() { echo "[$(date '+%F %T')] $*"; }
[ "$(git -C "$MM" rev-parse --short HEAD)" = 6409f5d ] || { echo "MicrobeMixer is not at 6409f5d"; exit 1; }
for pr in $PR; do
  prof=${pr%%.*}
  M=$OUT/mm.$prof
  if [ ! -d "$M" ]; then
    mkdir -p "$M"; cp -r "$MM/src" "$M/src"
    sed -i "s/art_illumina -ss MSv1 /art_illumina -ss $prof /" "$M/src/simulate.py"
    [ "$(grep -c "art_illumina -ss $prof " "$M/src/simulate.py")" = 1 ] || { echo "could not set the profile $prof"; exit 1; }
  fi
  for r in $REGIONS; do
    d=$OUT/$pr/$r; c=$d/credit; [ -s "$d/done" ] && continue
    rm -rf "$d"; mkdir -p "$c" "$SAMPLES/$pr"
    p=${r/V4_V4/V4}
    stamp "$pr $r: simulating"
    $CAP python3 "$M/src/sim_16s_reads.py" simulate --biome-abundance "${AB//%r/$r}" --primers "$MM/data/${p}_primers.txt" --silva-taxonomy "$TAX" \
        --silva-ref "$FASTA" --num-reads 10000000 --output-name "$r" --temp-dir "$d/" > "$d/stdout.log" 2>&1 \
      || { echo "FAILED: the simulation of $pr $r (see $d/stdout.log)"; exit 1; }
    rm -f "$d"/genus_*
    for m in 1 2; do
      awk 'int((NR-1)/4) % 50 == 0' "$d/${r}_mate_$m.fq" | gzip -1 > "$SAMPLES/$pr/${r}_$m.fq.gz"
      rm -f "$d/${r}_mate_$m.fq"
      stamp "$pr $r: mate $m"
      $CAP "$K" find --preset silva -j $THREADS "$IDX" "$SAMPLES/$pr/${r}_$m.fq.gz" "$d/$m.out" 2> "$d/find_$m.log" \
        || { echo "FAILED: find on $pr $r"; exit 1; }
    done
    for spec in lca:lca: one:one: cliffs:cliffs: list:list: norm:norm: freq:freq: prior:prior: mix:mix: em:em: \
                bayesflat:bayes:--prior+flat+--epsilon+1e-6,1e-4,1e-5,1e-7,1e-8 bayesauto:bayesauto: bayesem:bayesem:--em-iters+20; do
      name=${spec%%:*}; rest=${spec#*:}; mode=${rest%%:*}; opts=${rest#*:}; opts=${opts//+/ }
      stamp "$pr $r: $name"
      $CAP python3 "$HERE/scripts/score_silva.py" --mode $mode $opts --levels genus --nseq "$NSEQ" --mate1-listings "$d/1.out" --mate2-listings "$d/2.out" \
          --doc-id-to-traversal "$D2T" --silva-tax-ranks "$TAX" --readset-truthset "$d/${r}_seqtax.txt" --trav-to-length "$T2L" --exp-src "$EXP" \
          --output "$c/$name.$r.csv" --threads $THREADS > "$c/$name.$r.csv.log" 2>&1 || { echo "FAILED: scoring $name on $pr $r"; exit 1; }
    done
    rm -f "$d/1.out" "$d/2.out"; date > "$d/done"
  done
done
stamp "done: $PR"
