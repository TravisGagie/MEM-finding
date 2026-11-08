#!/usr/bin/env python3
"""Summarises the sweep's logs (movi7.sh): t/<config>.r<rep>/<file>.log.  Time per read is the whole processing time over the
200,000 reads of each of the 8 files, for both tools (Movi 2's own figure, and katka's total time), excluding loading."""
import glob, os, re, statistics as st, sys
D = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/research/katka/work/pc"); M = D + "/movi"
Z = lambda f: os.path.getsize(f) / 1e6
I = Z(M + "/idx/index.movi"); F = Z(M + "/idx/ftab.10.bin"); C = Z(D + "/bac.csa"); ROW = 58.518
print("Movi 2 index %.1f MB, ftab %.1f MB, row starts %.1f MB; RLCSA %.1f MB" % (I, F, ROW, C))
rows = [("ours", "RLCSA, full tags", D + "/bac.tag"), ("full", "Movi 2, full tags", D + "/bac.tag"), ("s2", "Movi 2, s = 2", D + "/bac.s2.tag"), ("s4", "Movi 2, s = 4", D + "/bac.s4.tag"),
        ("s8", "Movi 2, s = 8", D + "/bac.s8.tag"), ("s16", "Movi 2, s = 16", D + "/bac.s16.tag"), ("grammar", "Movi 2, grammar", D + "/g.tag"), ("bigrepair", "Movi 2, bigrepair", D + "/bg.tag")]
print("%-20s %8s %11s %8s %10s %24s %7s" % ("configuration", "tags MB", "boundaries", "labels", "index MB", "us per read, median (range)", "peak GB"))
for c, name, f in rows:
    t = []; pk = 0; b = l = None
    for d in sorted(glob.glob("%s/t/%s.r?" % (M, c))):
        w = 0.0
        for g in glob.glob(d + "/*.log"):
            s = open(g, errors="ignore").read()
            x = re.search(r"processing the reads: ([\d.]+) seconds", s) if c != "ours" else re.search(r"timing: .*\(total ([\d.]+) s\)", s)
            w += float(x[1]); p = re.search(r"peak (\d+) KB", s); pk = max(pk, int(p[1]) if p else 0)
            y = re.search(r"tags .*?: boundaries ([\d.]+) MB, labels ([\d.]+) MB", s) or re.search(r"tag array: .*?boundaries ([\d.]+) MB, labels ([\d.]+) MB", s)
            if y: b, l = y[1], y[2]
        t.append(w / 1.6)
    idx = Z(f) + (C if c == "ours" else I + F + ROW)
    print("%-20s %8.1f %11s %8s %10.1f %9.2f (%.2f-%.2f) %7.2f  n=%d" % (name, Z(f), b, l, idx, st.median(t), min(t), max(t), pk * 1024 / 1e9, len(t)))
