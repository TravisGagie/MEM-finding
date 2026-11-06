#!/usr/bin/env python3
"""Summarises experiments/narrow_timing.sh's logs (OUT) and experiments/csa_sampled_tags.sh's (CSATAGS).

    narrow_tables.py OUT CSATAGS IDX GRAMMAR
"""
import glob, os, re, statistics as st, sys
O, C, P, G = sys.argv[1:5]
def rep(pat):
    m = []; t = []; pk = 0; tab = None; tg = None
    for d in sorted(glob.glob(pat)):
        w = n = tt = 0
        for f in glob.glob(d + "/*.log"):
            s = open(f).read(); x = re.search(r"timing: (\d+) reads, 1 threads; mean per-read time ([\d.]+) us; total time / reads ([\d.]+) us", s)
            if not x: continue
            w += int(x[1]) * float(x[2]); tt += int(x[1]) * float(x[3]); n += int(x[1])
            p = re.search(r"peak (\d+) KB", s); pk = max(pk, int(p[1]) if p else 0)
            y = re.search(r"lookup table: (.*)", s); tab = y[1] if y else tab
            y = re.search(r"tag array: .*?boundaries ([\d.]+) MB, labels ([\d.]+) MB", s); tg = y.groups() if y else tg
        if n: m.append(w / n); t.append(tt / n)
    return m, t, pk, tab, tg
print("== the lookup table at every step: mean per read (us), median (range) of 5 repeats; peak GB")
names = [("cur", "current search, plain table"), ("c", "compact table"), ("ca", "compact, (a)"), ("cb", "compact, (b)"), ("cab", "compact, (a) and (b)"), ("pa", "plain, (a)"), ("pb", "plain, (b)"), ("pab", "plain, (a) and (b)")]
for k in (8, 10, 11, 12):
    print("k = %d" % k)
    for c, n in names:
        m, t, pk, tab, _ = rep("%s/%s%d.r?" % (O, c, k))
        if m: print("  %-30s %7.2f (%.2f-%.2f)  peak %.2f  n=%d  table: %s" % (n, st.median(m), min(m), max(m), pk * 1024 / 1e9, len(m), tab or ""))
print("\n== explicit RLCSA, no RLBWT loaded: tags")
Z = lambda f: os.path.getsize(f) / 1e6
for n, f in (("full", P + ".tag"), ("s2", P + ".s2.tag"), ("s4", P + ".s4.tag"), ("s8", P + ".s8.tag"), ("grammar", G + "/g.tag"), ("bigrepair", G + "/bg.tag")):
    m, t, pk, _, tg = rep("%s/%s.r?" % (C, n))
    print("  %-10s tags %7.1f MB (boundaries %s, labels %s)  index %7.1f MB  %7.2f (%.2f-%.2f) us  total/reads %.2f  peak %.2f GB  n=%d" % (n, Z(f), tg[0], tg[1], Z(P + ".csa") + Z(f), st.median(m), min(m), max(m), st.median(t), pk * 1024 / 1e9, len(m)))
