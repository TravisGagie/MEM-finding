#!/usr/bin/env python3
"""Naive Bayes with an EM prior beside the other schemes, on every read set (score_silva.py --mode bayesem).

    bayesem_tables.py [QUEUE] [TOOL_SILVA_DIR]

An accuracy is strict, at genus level, the mean of the four regions; an entry is the mean of the repeats or draws (range).
The dial's second block excludes reads from genera with ambiguous names; the other sets have none.
"""
import glob, os, re, sys
Q = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/tool/silva/queue"); A = os.path.expanduser(sys.argv[2] if len(sys.argv) > 2 else "~/tool/silva")
R = "V1_V2 V3_V4 V4_V4 V4_V5".split()
def load(ds):
    T = {}
    for d in ds:
        for f in glob.glob(d + "/*.csv"):
            m = re.search(r"(V\d_V\d)\.csv$", f)
            if not m or re.search(r"freq\.(grammar|bigrepair)", f): continue
            for l in open(f):
                p = l.strip().split(",")
                if len(p) < 4 or p[2] != "genus": continue
                k = dict(x.split("=") for x in p[4:] if "=" in x)
                if "strict_correct" in k: T.setdefault(p[1], {})[m[1]] = (int(k["strict_correct"]), int(k["reads"]), int(k["ambiguous_truth"]))
    return T
S = {"aquatic, MSv1": [load([Q + "/credit"])] + [load(glob.glob("%s/art_profiles/MSv1.%d/*/credit" % (A, k))) for k in (2, 3)],
     "aquatic, MSv3": [load(glob.glob("%s/art_profiles/MSv3.%d/*/credit" % (A, k))) for k in (1, 2, 3)]}
for a in ("1", "0", "-1"): S["dense dial, alpha = " + a] = [load(["%s/dial/a%s/credit" % (Q, a)])]
for s in ("uniform", "small"): S["sparse, " + s] = [load(glob.glob("%s/sparse/%s/MSv1.%d/*/credit" % (A, s, k))) for k in (1, 2, 3)]
K = [("freq", "proportional"), ("norm", "normalised"), ("em20", "EM 20"), ("bayesflat", "NB eps=1e-6"), ("bayesauto", "bayesauto"), ("bayesem5", "NB+EM 5"), ("bayesem20", "NB+EM 20"), ("bayesemmix5", "mix+EM 5"), ("bayesemmix20", "mix+EM 20")]
def m(T, k, ex):
    t = T.get(k, {})
    if not all(r in t for r in R): return None
    return sum(100 * ((t[r][0] - t[r][2]) / (t[r][1] - t[r][2]) if ex else t[r][0] / t[r][1]) for r in R) / 4
for ex, names in ((False, list(S)), (True, [n for n in S if "dense" in n])):
    print("== %s" % ("excluding reads from ambiguous genera" if ex else "all reads"))
    print("  %-24s" % "" + "".join("%-22s" % n for _, n in K))
    for s in names:
        row = ""
        for k, _ in K:
            v = [x for x in (m(T, k, ex) for T in S[s]) if x is not None]
            row += "%-22s" % (("%.2f" % (sum(v) / len(v)) + (" (%.2f-%.2f)" % (min(v), max(v)) if len(v) > 1 else "")) if v else "-")
        print("  %-24s" % s + row)
    print()
