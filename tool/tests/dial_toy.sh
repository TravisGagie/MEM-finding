#!/bin/bash
# The credit-scheme experiment (silva_queue.sh dial) end to end on a toy SILVA (tests/make_dial_toy.py): reads with
# k_g ~ N_g^alpha for alpha = 1, 0, -1, every scheme, the table, and the alpha fit.
#
#   tests/dial_toy.sh MICROBEMIXER EXP DOC_TO_TRAVERSAL WORK
#
# MICROBEMIXER: a clone of MicrobeMixer (with data/tax_slv_ssu_138.1.txt); EXP: cliffy-experiments/src;
# DOC_TO_TRAVERSAL: real SILVA genus paths (at least 12, e.g. katka/silva/work/ref/rna/doc_to_traversal.txt).
# Takes about a minute (2,400 read pairs per region and alpha).
set -euo pipefail
MM=$(cd "$1" && pwd); EXP=$(cd "$2" && pwd); D2T=$3; mkdir -p "$4"; W=$(cd "$4" && pwd)
HERE=$(cd "$(dirname "$0")/.." && pwd); K=$HERE/bin/katka
rm -rf "$W"/*
python3 "$HERE/tests/make_dial_toy.py" "$D2T" "$W"
(cd "$W" && $K prep list.tsv bac > prep.log 2>&1 && $K bwt bac > bwt.log 2>&1 && $K index --preset silva bac > index.log 2>&1)
S=$W/S; mkdir -p "$S/work/ref/rna" "$S/exp1_data"
echo "12 genera" > "$S/work/ref/done"; ln -sfn "$W/ref" "$S/work/ref/dna"
cp "$W/doc_to_traversal.txt" "$S/work/ref/rna/"; cp "$W/trav_to_length.txt" "$S/work/ref/"; cp "$MM/data/tax_slv_ssu_138.1.txt" "$S/exp1_data/"
Q=$W/q S=$S IDX=$W/bac EXP=$EXP MM=$MM FASTA=$W/silva.fasta NREADS=2400 THREADS=${THREADS:-4} PIN=0 \
  bash "$HERE/experiments/silva_queue.sh" dial > "$W/dial.log" 2>&1 || { tail "$W/dial.log"; exit 1; }
cmp -s "$W/q/nseq.txt" "$W/nseq.txt" || { echo "the runner's N_g table differs from the toy's"; exit 1; }
python3 "$HERE/scripts/dial_tables.py" "$W/q" | tee "$W/tables.txt" | sed -n '/mean over the regions/,/^$/p'
for a in 1 0 -1; do python3 "$HERE/scripts/silva_alpha_fit.py" "$W/doc_to_traversal.txt" "$W/nseq.txt" "$W/q/dial/a$a/V4_V4.seqtax.txt"; done
# what the toy must show: at alpha = -1 (small genera common) proportional credit loses to normalised credit, and EM
# recovers most of the difference; the endpoints of the mixtures equal their schemes
python3 - "$W/q" <<'PY'
import sys
Q = sys.argv[1]; regions = "V1_V2 V3_V4 V4_V4 V4_V5".split()
def acc(al, mode, label):
    v = []
    for r in regions:
        for l in open("%s/dial/a%s/credit/%s.%s.csv" % (Q, al, mode, r)):
            t = l.strip().split(","); kv = dict(x.split("=", 1) for x in t[4:])
            if t[1] == label and t[2] == "genus": v.append(int(kv["strict_correct"]) / int(kv["reads"]))
    return sum(v) / len(v)
fails = 0
def check(what, ok):
    global fails; print("%-70s %s" % (what, "ok" if ok else "FAIL")); fails += not ok
for al in ("1", "0", "-1"):
    check("alpha = %s: prior0 = mix0 = proportional" % al, acc(al, "prior", "prior0") == acc(al, "mix", "mix0") == acc(al, "freq", "freq"))
    check("alpha = %s: prior1 = normalised, mix1 = even" % al, acc(al, "prior", "prior1") == acc(al, "norm", "norm") and acc(al, "mix", "mix1") == acc(al, "list", "list"))
check("alpha = 0: naive Bayes with the flat prior = the reference (same prior)", acc("0", "bayesflat", "bayesflat") == acc("0", "bayes", "bayes"))
check("alpha = 1: naive Bayes with the database prior = the reference (same prior)", acc("1", "bayesdb", "bayesdb") == acc("1", "bayes", "bayes"))
check("alpha = -1: normalised beats proportional", acc("-1", "norm", "norm") > acc("-1", "freq", "freq"))
check("alpha = -1: EM (20) beats proportional", acc("-1", "em", "em20") > acc("-1", "freq", "freq"))
sys.exit(1 if fails else 0)
PY
echo "dial toy: ok ($W)"
