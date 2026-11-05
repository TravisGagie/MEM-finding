#!/usr/bin/env python3
"""Accuracy of every credit scheme on the aquatic reads by ART error profile (experiments/aquatic_art_profiles.sh).

    art_profile_tables.py [ART_DIR] [QUEUE]

MSv1.1 is the queue's existing reads (MicrobeMixer's own -ss MSv1), scored in QUEUE/credit.  Per repeat, an accuracy
is the mean of the four regions; each table shows the profile's mean over repeats and the range of the repeats.
strict: a tie for the top counts as wrong.  theirs: Cliffy's scoring, as in its classify_silva_readset.py.
"""
import glob, os, re, sys
A = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/tool/silva/art_profiles")
Q = os.path.expanduser(sys.argv[2] if len(sys.argv) > 2 else "~/tool/silva/queue")
REGIONS = "V1_V2 V3_V4 V4_V4 V4_V5".split()
NAMES = [("lca", "range (Cliffy's LCA query)"), ("one", "one label per MEM"), ("cliffs", "sample (simulated cliffs)"), ("list", "even"),
         ("norm", "normalised"), ("freq", "proportional"), ("prior0.5", "prior mixture, beta = 0.5"), ("mix0.5", "proportional/even mix, beta = 0.5"),
         ("em20", "EM, 20 iterations"), ("bayesflat", "naive Bayes, flat, eps = 1e-6"), ("bayesauto", "naive Bayes, mixture, eps at the peak")]
def load(d):
    T = {}
    for f in glob.glob(os.path.join(d, "*.csv")):
        m = re.search(r"(V\d_V\d)\.csv$", f)
        if not m or re.search(r"freq\.(grammar|bigrepair)", f): continue
        for l in open(f):
            p = l.strip().split(",")
            if len(p) < 4 or p[2] != "genus": continue
            kv = dict(x.split("=") for x in p[4:] if "=" in x)
            if "strict_correct" in kv:
                T.setdefault(p[1], {})[m.group(1)] = (int(kv["strict_correct"]) / int(kv["reads"]), float(kv["theirs"]))
    return T
R = {"MSv1.1": load(os.path.join(Q, "credit"))}
for pr in sorted(os.listdir(A)):
    if re.fullmatch(r"MSv\d\.\d", pr):
        t = {}
        for r in REGIONS:
            for k, v in load(os.path.join(A, pr, r, "credit")).items(): t.setdefault(k, {}).update(v)
        R[pr] = t
PROFS = sorted({p.split(".")[0] for p in R})
reps = {p: sorted(x for x in R if x.startswith(p + ".")) for p in PROFS}
print("repeats:", ", ".join("%s: %d" % (p, len(reps[p])) for p in PROFS), "(MSv1.1 is the existing reads)\n")
def val(pr, k, i, r=None):
    t = R[pr].get(k, {})
    if r: return 100 * t[r][i] if r in t else None
    return 100 * sum(t[x][i] for x in REGIONS) / 4 if all(x in t for x in REGIONS) else None
def cell(p, k, i, r=None):
    v = [val(x, k, i, r) for x in reps[p]]; v = [x for x in v if x is not None]
    if not v: return "-"
    return "%.2f (%.2f-%.2f; n=%d)" % (sum(v) / len(v), min(v), max(v), len(v))
for i, nm in ((0, "strict"), (1, "Cliffy's scoring")):
    print("== genus accuracy, %s, mean of the four regions: mean of the repeats (range)" % nm)
    print("  %-40s" % "scheme" + "".join("%-30s" % p for p in PROFS))
    for k, n in NAMES: print("  %-40s" % n + "".join("%-30s" % cell(p, k, i) for p in PROFS))
    print()
print("== range (Cliffy's LCA query) by region")
print("  %-16s" % "" + "".join("%-30s" % p for p in PROFS))
for i, nm in ((0, "strict"), (1, "Cliffy's")):
    for r in REGIONS: print("  %-6s %-9s" % (r, nm) + "".join("%-30s" % cell(p, "lca", i, r) for p in PROFS))
