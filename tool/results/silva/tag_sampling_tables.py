#!/usr/bin/env python3
"""Summarises experiments/tag_sampling.sh's logs: the median and range over the repeats of the mean per-read time, one thread.

    tag_sampling_tables.py OUT IDX
"""
import glob, os, re, statistics as st, sys
O, P = sys.argv[1], sys.argv[2]
Z = lambda f: os.path.getsize(f) / 1e6
rows = [("explicit RLCSA, full tags", "csa", [P + ".rlbwt", P + ".csa", P + ".tag"])]
for s in (1, 2, 4, 8):
    rows.append(("RLBWT, full tags" if s == 1 else "RLBWT, s = %d" % s, "s%d" % s, [P + ".rlbwt", P + (".tag" if s == 1 else ".s%d.tag" % s)]))
print("%-28s %10s %10s %24s %24s %9s" % ("configuration", "tags MB", "index MB", "mean per read, us", "total/reads, us", "peak GB"))
for name, c, fs in rows:
    m = []; t = []; pk = 0
    for d in sorted(glob.glob("%s/%s.r*" % (O, c))):
        w = n = tt = 0
        for f in glob.glob(d + "/*.log"):
            s = open(f).read(); x = re.search(r"timing: (\d+) reads, 1 threads; mean per-read time ([\d.]+) us; total time / reads ([\d.]+) us", s)
            w += int(x[1]) * float(x[2]); tt += int(x[1]) * float(x[3]); n += int(x[1])
            p = re.search(r"peak (\d+) KB", s); pk = max(pk, int(p[1]) if p else 0)
        m.append(w / n); t.append(tt / n)
    print("%-28s %10.1f %10.1f %9.2f (%.2f-%.2f) %5d %9.2f (%.2f-%.2f) %5.2f  n=%d" % (name, Z(fs[-1]), sum(Z(f) for f in fs), st.median(m), min(m), max(m), 0, st.median(t), min(t), max(t), pk * 1024 / 1e9, len(m)))
