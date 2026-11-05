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
# find --classify bayes gives the scorer's naive-Bayes calls (flat prior), read by read, over both mates and the scored
# documents.  The scorer's --exclude-out gives the documents; a toy from a real SILVA list has an "uncultured" genus
for al in 1 0 -1; do for e in 1e-6 1e-3; do
  A=$W/q/dial/a$al; r=V4_V4
  python3 "$HERE/scripts/score_silva.py" --mode bayes --prior flat --epsilon $e --levels genus --mate1-listings "$A/out/${r}_1.out" \
      --mate2-listings "$A/out/${r}_2.out" --doc-id-to-traversal "$S/work/ref/rna/doc_to_traversal.txt" --silva-tax-ranks "$S/exp1_data/tax_slv_ssu_138.1.txt" \
      --readset-truthset "$A/$r.seqtax.txt" --trav-to-length "$S/work/ref/trav_to_length.txt" --exp-src "$EXP" --nseq "$W/q/nseq.txt" \
      --output "$W/s.csv" --calls "$W/scorer_calls.tsv" --exclude-out "$W/excluded.txt" --threads ${THREADS:-4} > /dev/null 2>&1
  $K find --preset silva -j 4 --classify bayes --nseq "$W/q/nseq.txt" --epsilon $e --mate2 "$A/work/$r/${r}_mate_2.fq" --exclude-docs "$W/excluded.txt" \
      "$W/bac" "$A/work/$r/${r}_mate_1.fq" "$W/calls.tsv" 2> /dev/null
  cmp -s "$W/calls.tsv" "$W/scorer_calls.tsv" || { echo "FAIL: find --classify bayes differs from the scorer (alpha = $al, epsilon = $e)"; exit 1; }
  echo "find --classify bayes = the scorer's calls: alpha = $al, epsilon = $e, $(wc -l < "$W/calls.tsv") read pairs, $(wc -l < "$W/excluded.txt") excluded documents"
done; done
# naive Bayes with an EM prior against an independent implementation, on the toy's real listings, call for call: the
# scorer's calls must be among the independent top genera (within rounding), at each alpha
for al in 1 0 -1; do for me in "" "--mixture-eps 0.1"; do
  A=$W/q/dial/a$al; r=V4_V4
  python3 "$HERE/scripts/score_silva.py" --mode bayesem --em-iters 20 $me --levels genus --mate1-listings "$A/out/${r}_1.out" \
      --mate2-listings "$A/out/${r}_2.out" --doc-id-to-traversal "$S/work/ref/rna/doc_to_traversal.txt" --silva-tax-ranks "$S/exp1_data/tax_slv_ssu_138.1.txt" \
      --readset-truthset "$A/$r.seqtax.txt" --trav-to-length "$S/work/ref/trav_to_length.txt" --exp-src "$EXP" --nseq "$W/q/nseq.txt" \
      --output "$W/s.csv" --calls "$W/bem_calls.tsv" --exclude-out "$W/excluded.txt" --threads ${THREADS:-4} > /dev/null 2>&1 || { echo "FAIL: bayesem"; exit 1; }
  python3 - "$A/out/${r}_1.out" "$A/out/${r}_2.out" "$W/q/nseq.txt" "$W/excluded.txt" "$W/bem_calls.tsv" "${me##* }" <<'PY' || exit 1
import math, re, sys
m1, m2, nseqf, exf, callsf, mix = sys.argv[1:7]
mix = float(mix) if mix else None; eps = 1e-6
nseq = [int(x) for x in open(nseqf)]; G = len(nseq); N = sum(nseq); ex = {int(x) for x in open(exf)}
def rd(f):
    out = []
    for l in open(f):
        if l[0] == ">": out.append([l[1:].strip(), []])
        else:
            for mem in re.findall(r"\{([^}]*)\}", l):
                out[-1][1].append({int(x.split(":")[0]): int(x.split(":")[1]) for x in mem.split(",") if x})
    return out
A1, A2 = rd(m1), rd(m2); reads = []
for (n, f1), (_, f2) in zip(A1, A2):
    F = f1 + f2; S = {}
    for g in {g for f in F for g in f}:
        v = 0.0
        for f in F:
            if mix is None: v += math.log(f[g] / (2 * nseq[g])) if g in f else math.log(eps)
            else:
                q = mix * sum(f.values()) / (2 * N)
                v += math.log((1 - mix) * f.get(g, 0) / (2 * nseq[g]) + q) - math.log(q)
        S[g] = v
    reads.append((n, S))
pi = [1 / G] * G
for _ in range(20):
    T = [0.0] * G
    for _, S in reads:
        if not S: continue
        top = max(S[g] + math.log(pi[g]) for g in S if pi[g] > 0)
        w = {g: math.exp(S[g] + math.log(pi[g]) - top) for g in S if pi[g] > 0}; t = sum(w.values())
        for g, x in w.items(): T[g] += x / t
    tt = sum(T); pi = [x / tt for x in T]
bad = n = 0
for (name, S), l in zip(reads, open(callsf)):
    got = l.split()[1]; got = [] if got == "-" else [int(x) for x in got.split(",")]
    sc = {g: S[g] + math.log(pi[g]) for g in S if pi[g] > 0 and g not in ex}; n += 1
    top = max(sc.values()) if sc else 0
    near = {g for g, v in sc.items() if v > top - 1e-7}
    if not set(got) <= near or (near and not got): bad += 1
print("bayesem = an independent implementation, call for call: %d reads, %d disagree (%s)" % (n, bad, "mixture" if mix else "epsilon = 1e-6"))
sys.exit(1 if bad else 0)
PY
done; done
echo "dial toy: ok ($W)"
