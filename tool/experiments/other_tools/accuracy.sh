#!/bin/bash
# Genus accuracy of Kraken 2 or Tagger on the same read sets as KATKA's rows, 200,000 pairs per region, both mates.
#   accuracy.sh OUT kraken2 DB SET...         SET: MSv1.1 (the aquatic reads used everywhere), MSv1.2, MSv1.3, MSv3.1..3
#   accuracy.sh OUT tagger:L DIR SET...            (ART profile repeats), heldout1..3 (the held-out splits' reads)
# Scored with score_silva.py --mode list on listings that give each mate its tool's genus (to_listings.py): strict
# (a tie counts as wrong) and theirs (Cliffy's rule: ties, and reads without a genus, count as correct).  Tagger gives
# each mate a genus, so mates that disagree tie; Kraken 2 gives the pair one taxon, so it never ties, and Cliffy's
# rule would only count its unassigned pairs as correct.  Kraken 2 is also scored with Cliffy's own Kraken 2 scoring
# (classify_silva_readset.py's KrakenReadSet: with and without "vague positives").  Runs with THREADS (12) threads.
# Results: OUT/<SET>/<region>/{list.csv, kraken_cliffy.csv, *.log}.
set -uo pipefail
OUT=$1; TOOL=$2; IDX=$3; shift 3
HERE=$(cd "$(dirname "$0")" && pwd); TOP=$(cd "$HERE/../../.." && pwd); S=$TOP/katka/silva; EXP=$TOP/katka/cliffy-experiments/src
TH=${THREADS:-12}; X=$HOME/ext/moved_2026-11-05; TG=${TG:-$HOME/research/deps/tagger}; TAX=$S/exp1_data/tax_slv_ssu_138.1.txt
D2T=$S/work/ref/rna/doc_to_traversal.txt; T2L=$S/work/ref/trav_to_length.txt
CAP="systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0"
TMP=${TMPDIR:-/tmp}/other_tools.$$; mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
reads() {   # set region mate -> path
  case $1 in
    MSv1.1) echo $S/work/compare/reads_every50/${2}_$3.fq ;;
    MSv*) echo $X/art_profiles_samples/$1/${2}_$3.fq.gz ;;
    heldout1) echo $X/heldout_samples/MSv1.1/${2}_$3.fq.gz ;;
    heldout*) echo $X/heldout_samples${1#heldout}/MSv1.1/${2}_$3.fq.gz ;;
  esac; }
truth() {
  case $1 in
    MSv1.1) echo $S/work/reads/aquatic/$2/${2}_seqtax.txt ;;
    MSv*) echo $HOME/research/katka/work/art_profiles/$1/$2/${2}_seqtax.txt ;;
    heldout1) echo $HOME/research/katka/work/heldout/runs/MSv1.1/$2/${2}_seqtax.txt ;;
    heldout*) echo $HOME/research/katka/work/heldout${1#heldout}/runs/MSv1.1/$2/${2}_seqtax.txt ;;
  esac; }
for set in "$@"; do for r in V1_V2 V3_V4 V4_V4 V4_V5; do
  d=$OUT/$set/$r; [ -s $d/list.csv ] && continue; mkdir -p $d
  for m in 1 2; do f=$(reads $set $r $m); [ -s "$f" ] || { echo "no reads $f"; exit 1; }
    case $f in *.gz) zcat $f > $TMP/$m.fq ;; *) cp $f $TMP/$m.fq ;; esac; done
  T=$(truth $set $r); [ -s "$T" ] || { echo "no truth $T"; exit 1; }
  case $TOOL in
    kraken2)
      $CAP kraken2 --db $IDX --threads $TH --paired $TMP/1.fq $TMP/2.fq > $TMP/k.out 2> $d/kraken2.log || { echo "kraken2 failed $set $r"; exit 1; }
      python3 $HERE/to_listings.py kraken2 $TMP/k.out $TMP/1.fq $TMP/l1 $TMP/l2 $TAX $D2T 2> $d/convert.log
      mkdir -p $d/cliffy; python3 - "$EXP" "$TAX" "$T" "$TMP/k.out" "$d/cliffy/" > $d/cliffy.log 2>&1 <<'PY'
import sys, importlib.util
exp, tax, truth, kout, outdir = sys.argv[1:6]
spec = importlib.util.spec_from_file_location("csr", exp + "/classify_silva_readset.py")
csr = importlib.util.module_from_spec(spec); spec.loader.exec_module(csr)
csr.KrakenReadSet(tax, truth, kout, None).generate_sensitivity_plot(outdir, "output")
PY
      grep ',genus,' $d/cliffy/output.classification_results.csv > $d/kraken_cliffy.csv; rm -f $d/cliffy/TP_FP_VP_results.txt ;;
    tagger:*)
      L=${TOOL#tagger:}
      for m in 1 2; do $CAP $TG/build/tagger -f $TMP/$m.fq -r $IDX/silva -L $L -t $TH -R -o $TMP/t$m.out > $d/tagger_$m.log 2>&1 || { echo "tagger failed $set $r"; exit 1; }; done
      python3 $HERE/to_listings.py tagger $TMP/t1.out $TMP/t2.out $TMP/1.fq $TMP/2.fq $TMP/l1 $TMP/l2 ;;
  esac
  $CAP python3 $TOP/tool/scripts/score_silva.py --mode list --levels genus --mate1-listings $TMP/l1 --mate2-listings $TMP/l2 \
     --doc-id-to-traversal $D2T --silva-tax-ranks $TAX --readset-truthset $T --trav-to-length $T2L --exp-src $EXP \
     --output $d/list.csv --threads $TH > $d/score.log 2>&1 || { echo "scoring failed $set $r"; exit 1; }
  echo "$TOOL $set $r $(cut -d, -f4,10,11 $d/list.csv) $(cat $d/kraken_cliffy.csv 2>/dev/null | tr '\n' ' ')"
done; done
