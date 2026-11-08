#!/usr/bin/env python3
"""Naive Bayes with an EM prior learned from a sample of the reads, then every read classified with it (emsample.sh's CSVs).

    em_sample_table.py DIR

An entry is strict genus accuracy, the mean of the four regions, for each set; aquatic is one read set, uniform and small the mean of 3
draws; the range is over the 3 sampling seeds (for uniform and small, the seed's mean over the draws).  Also the loss against the full-set prior.
"""
import glob, os, re, sys, statistics as st
D = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/research/katka/work/emsample")
R = "V1_V2 V3_V4 V4_V4 V4_V5".split(); T = {}
for f in glob.glob(D + "/*.csv"):
    s, k, r = os.path.basename(f)[:-4].split(".")
    for l in open(f):
        p = l.strip().split(",")
        if len(p) > 3 and p[2] == "genus":
            kv = dict(x.split("=") for x in p[4:] if "=" in x)
            T.setdefault(s, {}).setdefault(k, {}).setdefault(p[1], {})[r] = (int(kv["strict_correct"]) / int(kv["reads"]) * 100, int(kv.get("unclassified", 0)))
def m(s, lab):
    v = []
    for k in T[s]:
        t = T[s][k].get(lab, {})
        if len(t) == 4: v.append(sum(t[r][0] for r in R) / 4)
    return sum(v) / len(v) if v else None
sets = [x for x in ("aquatic", "uniform", "small") if x in T]
print("full-set prior: " + "  ".join("%s %.2f" % (s, m(s, "bayesem20")) for s in sets) + "; flat prior, for reference, is in the other tables\n")
print("%-8s %-9s" % ("reads", "smoothing") + "".join("%-30s" % s for s in sets))
for n in (1000, 10000, 100000):
    for sm, name in (("0", "none"), ("p0.5", "0.5 reads"), ("p1", "1 read"), ("m0.01", "flat 1%"), ("m0.1", "flat 10%")):
        row = "%-8d %-9s" % (n, name)
        for s in sets:
            v = [m(s, "bayesem20n%ds%d%s" % (n, q, sm)) for q in (1, 2, 3)]; f = m(s, "bayesem20")
            row += "%-30s" % ("%.2f (%.2f-%.2f) %+.2f" % (sum(v) / 3, min(v), max(v), sum(v) / 3 - f))
        print(row)
