#!/usr/bin/env python3
"""The sparse dial: every credit scheme on the aquatic abundances moved to random genera, beside the aquatic reads.

    sparse_tables.py [SPARSE_DIR] [ART_DIR] [QUEUE]

Per draw an accuracy is the mean of the four regions; a table entry is the mean of the draws (range).  aquatic is the
3 MSv1 repeats of art_profile_tables.py.  Strict genus accuracy; no read here is from a genus with an ambiguous name
(checked), so excluding them changes nothing.
"""
import glob, os, re, sys
A = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/tool/silva/sparse")
ART = os.path.expanduser(sys.argv[2] if len(sys.argv) > 2 else "~/tool/silva/art_profiles")
Q = os.path.expanduser(sys.argv[3] if len(sys.argv) > 3 else "~/tool/silva/queue")
REGIONS = "V1_V2 V3_V4 V4_V4 V4_V5".split()
NAMES = [("lca", "range (Cliffy's LCA query)"), ("one", "one label per MEM"), ("cliffs", "sample (simulated cliffs)"), ("list", "even"),
         ("norm", "normalised"), ("freq", "proportional"), ("prior0.1", "prior mixture, beta = 0.1"), ("prior0.5", "prior mixture, beta = 0.5"),
         ("prior0.75", "prior mixture, beta = 0.75"), ("mix0.5", "proportional/even mix, beta = 0.5"), ("em5", "EM, 5 iterations"), ("em20", "EM, 20 iterations"),
         ("bayesflat", "naive Bayes, flat, eps = 1e-6"), ("bayesauto", "naive Bayes, mixture, eps at the peak")]
amb = 0
def load(ds):
    global amb
    T = {}
    for d in ds:
        for f in glob.glob(os.path.join(d, "*.csv")):
            m = re.search(r"(V\d_V\d)\.csv$", f)
            if not m or re.search(r"freq\.(grammar|bigrepair)", f): continue
            for l in open(f):
                p = l.strip().split(",")
                if len(p) < 4 or p[2] != "genus": continue
                kv = dict(x.split("=") for x in p[4:] if "=" in x)
                if "strict_correct" in kv:
                    T.setdefault(p[1], {})[m.group(1)] = int(kv["strict_correct"]) / int(kv["reads"]); amb += int(kv["ambiguous_truth"])
    return T
SETS = {}
for s in ("uniform", "small"):
    SETS[s] = {k: load(glob.glob(os.path.join(A, s, k, "*", "credit"))) for k in sorted(os.listdir(os.path.join(A, s))) if re.fullmatch(r"MSv1\.\d", k)}
print("reads from ambiguous genera in the sparse draws: %d" % amb)
SETS["aquatic"] = {"MSv1.1": load([os.path.join(Q, "credit")])}
for k in ("MSv1.2", "MSv1.3"): SETS["aquatic"][k] = load(glob.glob(os.path.join(ART, k, "*", "credit")))
def mean(T, k, r=None):
    t = T.get(k, {})
    if r: return 100 * t[r] if r in t else None
    return 100 * sum(t[x] for x in REGIONS) / 4 if all(x in t for x in REGIONS) else None
def cell(s, k, r=None):
    v = [mean(T, k, r) for T in SETS[s].values()]; v = [x for x in v if x is not None]
    return "%6.2f (%.2f-%.2f; n=%d)" % (sum(v) / len(v), min(v), max(v), len(v)) if v else "-"
cols = ["aquatic", "uniform", "small"]
print("\n== strict genus accuracy, mean of the four regions: mean of the draws (range)")
print("  %-40s" % "scheme" + "".join("%-30s" % c for c in cols))
for k, n in NAMES: print("  %-40s" % n + "".join("%-30s" % cell(c, k) for c in cols))
for k, n in (("lca", "range"), ("freq", "proportional"), ("norm", "normalised"), ("em20", "EM, 20 iterations"), ("bayesflat", "naive Bayes, eps = 1e-6")):
    print("\n== %s, by region" % n); print("  %-8s" % "" + "".join("%-30s" % c for c in cols))
    for r in REGIONS: print("  %-8s" % r + "".join("%-30s" % cell(c, k, r) for c in cols))
print("\n== bayesauto's epsilon")
for c in ("uniform", "small"):
    for k in sorted(SETS[c]):
        e = []
        for r in REGIONS:
            f = os.path.join(A, c, k, r, "credit", "bayesauto.%s.csv.log" % r)
            m = re.findall(r"bayesauto: epsilon = ([\d.e-]+)", open(f).read()) if os.path.exists(f) else []
            e.append(m[-1] if m else "-")
        print("  %-8s %s: %s" % (c, k, " ".join(e)))
