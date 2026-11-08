#!/usr/bin/env python3
"""Every SILVA scheme, all reads and excluding reads from genera with ambiguous names; the scoring CSVs only.

    rescore_tables.py [Q]        Q: the queue directory (default ~/research/katka/work/queue)

A read whose true genus name is missing from SILVA's genus rank, or sits at more than one node of it, counts as
correct in the scorer (as in Cliffy's) and is counted in ambiguous_truth.  So, per region:
    all reads:  strict_correct / reads
    excluding:  (strict_correct - ambiguous_truth) / (reads - ambiguous_truth)
A table entry is the mean of the four regions.
"""
import glob, math, os, re, sys
Q = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/research/katka/work/queue")
REGIONS = "V1_V2 V3_V4 V4_V4 V4_V5".split()
NAMES = [("lca", "range (Cliffy's LCA query)"), ("one", "one label per MEM"), ("cliffs", "sample (simulated cliffs)"),
         ("list", "even"), ("norm", "normalised"), ("freq", "proportional"),
         ("prior0.1", "prior mixture, beta = 0.1"), ("prior0.25", "prior mixture, beta = 0.25"),
         ("prior0.5", "prior mixture, beta = 0.5"), ("prior0.75", "prior mixture, beta = 0.75"),
         ("mix0.1", "proportional/even mix, beta = 0.1"), ("mix0.25", "proportional/even mix, beta = 0.25"),
         ("mix0.5", "proportional/even mix, beta = 0.5"), ("mix0.75", "proportional/even mix, beta = 0.75"),
         ("em5", "EM, 5 iterations"), ("em20", "EM, 20 iterations"),
         ("bayes", "naive Bayes, the sample's prior, eps = 1e-6"), ("bayesdb", "naive Bayes, database prior, eps = 1e-6"),
         ("bayesflat0.01", "naive Bayes, flat, eps = 1e-2"), ("bayesflat0.001", "naive Bayes, flat, eps = 1e-3"),
         ("bayesflat0.0001", "naive Bayes, flat, eps = 1e-4"), ("bayesflat", "naive Bayes, flat, eps = 1e-6"),
         ("bayesflat1e-09", "naive Bayes, flat, eps = 1e-9"), ("bayesflat1e-12", "naive Bayes, flat, eps = 1e-12"),
         ("bayesflat1e-300", "naive Bayes, flat, eps = 1e-300 (check)"), ("bayesflat0", "naive Bayes, flat, the limit eps -> 0"),
         ("bayesauto", "naive Bayes, mixture, eps at the likelihood's peak")]
def load(d, skip=()):
    T = {}
    for f in glob.glob(os.path.join(d, "*.csv")):
        m = re.fullmatch(r".*?\.?(V\d_V\d)\.csv", os.path.basename(f))
        if not m or os.path.basename(f).split(".")[0] in skip or re.match(r"freq\.(grammar|bigrepair)", os.path.basename(f)): continue
        for l in open(f):
            p = l.strip().split(",")
            if len(p) < 3 or p[2] != "genus": continue
            kv = dict(x.split("=") for x in p[4:] if "=" in x)
            T.setdefault(p[1], {})[m.group(1)] = (int(kv["strict_correct"]), int(kv["reads"]), int(kv["ambiguous_truth"]))
    return T
def acc(T, k, excl):
    if k not in T or any(r not in T[k] for r in REGIONS): return None
    v = []
    for r in REGIONS:
        s, n, a = T[k][r]
        v.append(100 * ((s - a) / (n - a) if excl else s / n))
    return sum(v) / len(v)
def table(title, cols, Ts):
    print("== " + title)
    print("  %-52s" % "scheme" + "".join("%20s" % c for c in cols))
    print("  %-52s" % "" + "".join("%10s%10s" % ("all", "excl.") for c in cols))
    for k, name in NAMES:
        row = [(acc(T, k, False), acc(T, k, True)) for T in Ts]
        if all(x[0] is None for x in row): continue
        print("  %-52s" % name + "".join("%10s%10s" % tuple("-" if y is None else "%.2f" % y for y in x) for x in row))
    print()
def amb(T):
    t = T["freq"]; return [t[r][2] for r in REGIONS], t[REGIONS[0]][1]
D = [load(os.path.join(Q, "dial/a%s/credit" % a)) for a in ("1", "0", "-1")]
U = load(os.path.join(Q, "uniform/credit2")); A = load(os.path.join(Q, "credit"))
print("Strict genus accuracy (%), mean of the four regions: all reads, and excluding reads from genera with ambiguous names.\n")
print("Reads from genera with ambiguous names, of 200,000 pairs, by region (%s):" % " ".join(REGIONS))
for name, T in zip(("dial, alpha = 1", "dial, alpha = 0", "dial, alpha = -1", "uniform", "aquatic"), D + [U, A]):
    print("  %-18s %s" % (name, "  ".join("%6d" % x for x in amb(T)[0])))
print()
table("the dial", ["alpha = 1", "alpha = 0", "alpha = -1"], D)
table("the uniform simulation and the real aquatic reads", ["uniform", "aquatic"], [U, A])
# bayesauto's likelihood
print("== bayesauto: the log-likelihood's peak over the grid (in units of 10^6, at the peak and its neighbours)")
SETS = [("alpha = %s" % al, "dial/a%s/credit" % al) for al in ("1", "0", "-1")] + [("uniform", "uniform/credit2"), ("aquatic", "credit")]
for al, d in SETS:
    for r in REGIONS:
        s = open(os.path.join(Q, d, "bayesauto.%s.csv.log" % r)).read()
        m = re.search(r"log-likelihood by epsilon: (.*)", s); pts = [x.split(": ") for x in m.group(1).split(", ")]
        e = [float(x[0]) for x in pts]; ll = [float(x[1]) for x in pts]; k = max(range(len(ll)), key=lambda i: ll[i])
        near = "  ".join("%g: %.3f" % (e[i], ll[i] / 1e6) for i in range(max(0, k - 1), min(len(e), k + 2)))
        print("  %-11s %s  peak at eps = %-5g %s   (%s)" % (al, r, e[k], "AT THE EDGE" if k in (0, len(e) - 1) else "", near))
