#!/usr/bin/env python3
"""Tests of the credit schemes added to scripts/score_silva.py (prior and linear mixtures, EM) on hand-made listings.

    test_score_schemes.py EXP TAX DOC_TO_TRAVERSAL TRAV_TO_LENGTH [WORK]

EXP: cliffy-experiments/src; TAX: tax_slv_ssu_138.1.txt; DOC_TO_TRAVERSAL and TRAV_TO_LENGTH: a reference's files
with at least 4 genera (e.g. the toy's).  Three checks:
  1. endpoints, on 400 random reads with many exact ties: prior0 = freq, mix0 = freq, em0 = freq, prior1 = norm, mix1 = list, line for
     line at every level (all counts, with the label set aside);
  2. every scheme, genus level, against an independent implementation of its formula (exact rationals for prior,
     mix, norm, freq and list; floats for EM), on 400 other random reads whose genus totals don't tie: correct,
     unclassified and strictly correct counts must agree;
  3. hand-worked cases: one match {A:1, B:3} with N_A = 1, N_B = 9 and 12 genera (N = 20), truth A.  Proportional
     credit picks B; normalised credit picks A; the prior mixture picks A only for beta > 9/14 (0.75 and 1); the
     linear mixture never picks A alone (a tie at beta = 1, as even credit).  And a read {C:1, D:1} with
     N_C = N_D = 1, truth C, among 10 reads matching only C: proportional credit ties, EM learns that C is common
     and picks it.
"""
import os, random, subprocess, sys, tempfile
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
SCORE = os.path.join(HERE, "..", "scripts", "score_silva.py")
BETAS = [0, 0.1, 0.25, 0.5, 0.75, 1]
EXP, TAX, D2T, T2L = sys.argv[1:5]
W = sys.argv[5] if len(sys.argv) > 5 else tempfile.mkdtemp()
os.makedirs(W, exist_ok=True)
genera = [l.split(None, 1)[1].strip() for l in open(D2T) if l.strip()]
G = len(genera)
assert G >= 4, "need at least 4 genera"


def write_set(name, reads, nseq):
    """reads: (truth doc, mate-1 features, mate-2 features), a feature being (length, {doc: count})"""
    d = os.path.join(W, name); os.makedirs(d, exist_ok=True)
    with open(d + "/seqtax.txt", "w") as o:
        for i, g in enumerate(genera):
            o.write("%s\t%d\t%d\n" % (g, i + 1, sum(1 for r in reads if r[0] == i)))
    for m in (1, 2):
        with open("%s/m%d" % (d, m), "w") as o:
            for k, r in enumerate(reads):
                o.write(">genus_%d_read_%d_mate%d\n" % (r[0] + 1, k + 1, m))
                o.write(" ".join("[0,%d] {%s}" % (L - 1, ",".join("%d:%d" % (x, c) for x, c in sorted(f.items())))
                                 for L, f in r[m]) + "\n")
    open(d + "/nseq.txt", "w").write("".join("%d\n" % n for n in nseq))
    return d


def score(d, mode, extra=()):
    out = "%s/%s.csv" % (d, mode)
    cmd = [sys.executable, SCORE, "--mode", mode, "--mate1-listings", d + "/m1", "--mate2-listings", d + "/m2",
           "--doc-id-to-traversal", D2T, "--silva-tax-ranks", TAX, "--readset-truthset", d + "/seqtax.txt",
           "--trav-to-length", T2L, "--exp-src", EXP, "--output", out, "--nseq", d + "/nseq.txt", "--threads", "2"] + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit("score_silva.py --mode %s failed:\n%s" % (mode, r.stderr[-2000:]))
    lines = {}
    for l in open(out):
        t = l.strip().split(",")
        lines.setdefault(t[1], {})[t[2]] = t[3:]
    return lines


def reference(reads, nseq, scheme, beta=None, iters=None):
    """genus-level (correct, unclassified, strict) under the scheme, independently of score_silva.py"""
    N, Gn = sum(nseq), len(nseq)
    def weights(f, L, pi=None):
        docs = sorted(f)
        if scheme == "freq": x = {g: F(f[g]) for g in docs}
        elif scheme == "norm": x = {g: F(f[g], nseq[g]) for g in docs}
        elif scheme == "list": x = {g: F(1) for g in docs}
        elif scheme == "prior": x = {g: F(f[g], nseq[g]) * ((1 - F(beta)) * F(nseq[g], N) + F(beta) / Gn) for g in docs}
        elif scheme == "mix":
            tot = sum(f.values())
            return {g: L * ((1 - F(beta)) * F(f[g], tot) + F(beta) / len(docs)) for g in docs}
        elif scheme == "em": x = {g: f[g] * pi[g] / nseq[g] for g in docs}
        tot = sum(x.values())
        return {g: L * x[g] / tot for g in docs} if tot > 0 else {}
    pi = None
    if scheme == "em":
        pi = [n / N for n in nseq]
        for _ in range(iters):
            T = [0.0] * Gn
            for r in reads:
                for m in (1, 2):
                    for L, f in r[m]:
                        for g, w in weights(f, L, pi).items(): T[g] += w
            tt = sum(T); pi = [t / tt for t in T]
    tp = uncl = strict = 0
    for r in reads:
        acc = {}
        for m in (1, 2):
            for L, f in r[m]:
                for g, w in weights(f, L, pi).items(): acc[g] = acc.get(g, 0) + w
        if not acc: uncl += 1; continue
        top = max(acc.values())
        if acc.get(r[0]) == top:
            tp += 1; strict += sum(1 for v in acc.values() if v == top) == 1
    return tp, uncl, strict


def counts(fields):
    kv = dict(x.split("=", 1) for x in fields[1:])
    return int(kv["correct"]), int(kv["unclassified"]), int(kv["strict_correct"])


fails = 0
def check(what, ok):
    global fails
    print("%-74s %s" % (what, "ok" if ok else "FAIL")); fails += not ok


# 1. endpoints on random reads
rnd = random.Random(7)
nseq = [rnd.choice([1, 2, 3, 5, 9, 40]) for _ in range(G)]
def feat():
    docs = rnd.sample(range(G), rnd.choice([1, 1, 2, 2, 3, 4]))
    return (rnd.choice([30, 47, 60, 100]), {g: rnd.choice([1, 1, 2, 3, 7]) for g in docs})
reads = [(rnd.randrange(G), [feat() for _ in range(rnd.randrange(0, 4))], [feat() for _ in range(rnd.randrange(0, 4))])
         for _ in range(400)]
d = write_set("random", reads, nseq)
S = {m: score(d, m) for m in ("freq", "norm", "list")}
S.update(score(d, "prior")); S.update(score(d, "mix")); S.update(score(d, "em", ["--em-iters", "0,5,20"]))
for a, b in (("prior0", "freq"), ("mix0", "freq"), ("em0", "freq"), ("prior1", "norm"), ("mix1", "list")):
    B = S[b][b] if b in S[b] else S[b]
    check("endpoint: %s = %s at every level" % (a, b), S[a] == B)
# 2. against the reference implementation, on reads whose genus totals don't tie (float and exact sums then agree)
nseq = [rnd.randint(1, 500) for _ in range(G)]
def feat():
    docs = rnd.sample(range(G), rnd.choice([1, 2, 3, 4]))
    return (rnd.randint(30, 250), {g: rnd.randint(1, 1000) for g in docs})
reads = [(rnd.randrange(G), [feat() for _ in range(rnd.randrange(1, 4))], [feat() for _ in range(rnd.randrange(1, 4))])
         for _ in range(400)]
d = write_set("generic", reads, nseq)
BAYES_READS, BAYES_NSEQ = reads, nseq
S = {m: score(d, m) for m in ("freq", "norm", "list")}
S.update(score(d, "prior")); S.update(score(d, "mix")); S.update(score(d, "em", ["--em-iters", "0,5,20"]))
for m in ("freq", "norm", "list"):
    check("reference: %s" % m, counts(S[m][m]["genus"]) == reference(reads, nseq, m))
for b in BETAS:
    for m in ("prior", "mix"):
        check("reference: %s, beta = %g" % (m, b), counts(S["%s%g" % (m, b)]["genus"]) == reference(reads, nseq, m, beta=b))
for k in (0, 5, 20):
    check("reference: em, %d iterations" % k, counts(S["em%d" % k]["genus"]) == reference(reads, nseq, "em", iters=k))
# 3. hand-worked cases
nseq = [1, 9] + [1] * (G - 2)
reads = [(0, [(100, {0: 1, 1: 3})], [(100, {0: 1, 1: 3})])]
d = write_set("hand1", reads, nseq)
H = {m: score(d, m)[m]["genus"] for m in ("freq", "norm", "list")}
for m in ("prior", "mix"):
    H.update({k: v["genus"] for k, v in score(d, m).items()})
strict = lambda k: counts(H[k])[2]
correct = lambda k: counts(H[k])[0]
check("hand: proportional credit picks B (wrong)", strict("freq") == 0 and correct("freq") == 0)
check("hand: normalised credit picks A", strict("norm") == 1)
check("hand: prior mixture picks A for beta = 0.75 and 1 only (threshold 9/14)",
      [strict("prior%g" % b) for b in BETAS] == [0, 0, 0, 0, 1, 1])
check("hand: linear mixture never picks A alone; ties at beta = 1",
      [strict("mix%g" % b) for b in BETAS] == [0] * 6 and [correct("mix%g" % b) for b in BETAS] == [0] * 5 + [1])
check("hand: even credit ties (correct only with ties counted)", strict("list") == 0 and correct("list") == 1)
nseq = [1] * G
reads = [(2, [(80, {2: 1})], [(80, {2: 1})]) for _ in range(10)] + [(2, [(100, {2: 1, 3: 1})], [(100, {2: 1, 3: 1})])]
d = write_set("hand2", reads, nseq)
E = {k: v["genus"] for k, v in score(d, "em", ["--em-iters", "0,5,20"]).items()}
E["freq"] = score(d, "freq")["freq"]["genus"]
check("hand: proportional credit ties on {C, D} (10 of 11 strictly correct)", counts(E["freq"])[2] == 10)
check("hand: EM (5 and 20 iterations) learns C is common: 11 of 11", counts(E["em5"])[2] == 11 and counts(E["em20"])[2] == 11)
# the Bayes reference: a read with matches {A:1, B:5} and {A:1}, N_A = 1, N_B = 10.  Summed (proportional) credit
# gives B 5/6 of the first match's length (200) and A 1/6 + all of the second (100): B wins, 167 to 133.  The
# product gives A 1/2 * 1/2 (both strands), B 5/20 * epsilon: A wins, with any prior N_g^alpha, alpha in -1..1
nseq = [1, 10] + [1] * (G - 2)
reads = [(0, [(200, {0: 1, 1: 5})], [(100, {0: 1})])]
d = write_set("hand3", reads, nseq)
check("hand: proportional credit picks B for {A:1, B:5} + {A:1}", counts(score(d, "freq")["freq"]["genus"])[2] == 0)
check("hand: the Bayes reference picks A (alpha = 1, 0, -1)",
      all(counts(score(d, "bayes", ["--alpha", al])["bayes"]["genus"])[2] == 1 for al in ("1", "0", "-1")))
# naive Bayes with practical priors: flat = the sample prior with alpha = 0, db = with alpha = 1, line for line, on
# the generic reads; and against an independent implementation, for each epsilon
import math
EPS = ["1e-2", "1e-3", "1e-4", "1e-6", "1e-9"]
d = os.path.join(W, "generic")
gen_reads, gen_nseq = BAYES_READS, BAYES_NSEQ
flat = score(d, "bayes", ["--prior", "flat", "--epsilon", ",".join(EPS)])
dbp = score(d, "bayes", ["--prior", "db"])
s0 = score(d, "bayes", ["--prior", "sample", "--alpha", "0"])["bayes"]
s1 = score(d, "bayes", ["--prior", "sample", "--alpha", "1"])["bayes"]
check("bayes: the flat prior = the sample prior with alpha = 0 (epsilon 1e-6)", flat["bayesflat"] == s0)
check("bayes: the database prior = the sample prior with alpha = 1", dbp["bayesdb"] == s1)
check("bayes: one line per epsilon", sorted(flat) == sorted(["bayesflat0.01", "bayesflat0.001", "bayesflat0.0001", "bayesflat", "bayesflat1e-09"]))
def bayes_ref(reads, nseq, alpha, eps):
    tp = uncl = strict = 0
    for r in reads:
        feats = [f for m in (1, 2) for f in r[m]]
        cand = {g for _, f in feats for g in f}
        if not cand: uncl += 1; continue
        sc = {g: alpha * math.log(nseq[g]) + sum(math.log(f[g] / (2 * nseq[g])) if g in f else math.log(eps) for _, f in feats) for g in cand}
        top = max(sc.values())
        if sc.get(r[0]) == top:
            tp += 1; strict += sum(1 for v in sc.values() if v == top) == 1
    return tp, uncl, strict
for e, lab in zip(EPS, ["bayesflat0.01", "bayesflat0.001", "bayesflat0.0001", "bayesflat", "bayesflat1e-09"]):
    check("reference: bayes, flat prior, epsilon = %s" % e, counts(flat[lab]["genus"]) == bayes_ref(gen_reads, gen_nseq, 0, float(e)))
check("reference: bayes, database prior", counts(dbp["bayesdb"]["genus"]) == bayes_ref(gen_reads, gen_nseq, 1, 1e-6))
# epsilon decides: truth A (N_A = 1) has p = 1 on two matches and lacks a third; B (N_B = 10) has p = 1/20 on all three:
# A scores log epsilon, B 3 log(1/20) = log 1.25e-4, so A wins for epsilon = 1e-2 and 1e-3 only (flat prior)
nseq = [1, 10] + [1] * (G - 2)
reads = [(0, [(60, {0: 2, 1: 1}), (60, {0: 2, 1: 1})], [(60, {1: 1})])]
d = write_set("hand4", reads, nseq)
F4 = score(d, "bayes", ["--prior", "flat", "--epsilon", ",".join(EPS)])
check("hand: the flat-prior Bayes picks A for epsilon = 1e-2 and 1e-3 only",
      [counts(F4[l]["genus"])[2] for l in ["bayesflat0.01", "bayesflat0.001", "bayesflat0.0001", "bayesflat", "bayesflat1e-09"]] == [1, 1, 0, 0, 0])
print("%d failures" % fails)
sys.exit(1 if fails else 0)
