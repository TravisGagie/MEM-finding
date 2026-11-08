#!/usr/bin/env python3
"""Summarises movi13.sh's logs (t6/<config>.r<rep>/<file>.log): µs per read over the 8 files, one thread, setup excluded.
RLCSA index sizes include the 17.8 MB lookup table (-F 10, built at load); Movi's include its ftab (k = 10) and the row list."""
import glob, os, re, statistics as st, sys
D = os.path.expanduser("~/research/katka/work/pc"); M = D + "/movi"; T = sys.argv[1] if len(sys.argv) > 1 else "t6"
Z = lambda f: os.path.getsize(f) / 1e6
LT = 17.8
base = {"csa": Z(D + "/bac.csa") + LT, "ef": Z(D + "/bac.efcsa") + LT, "movi": Z(M + "/idx/index.movi") + Z(M + "/idx/ftab.10.bin") + 58.518}
tags = {"full": D + "/bac.tag", "grammar": D + "/g.tag", "s2": D + "/bac.s2.tag", "s4": D + "/bac.s4.tag", "s8": D + "/bac.s8.tag"}
names = {"csa": "explicit RLCSA", "ef": "Elias-Fano RLCSA", "movi": "Movi 2, row list"}
rows = []
for k, ts in (("csa", ("full", "grammar", "s2", "s4", "s8")), ("ef", ("full", "grammar", "s4", "s8")), ("movi", ("full", "grammar"))):
    for t in ts:
        c = "%s_%s" % (k, t); v = []; b = l = None
        for d in sorted(glob.glob("%s/%s/%s.r?" % (M, T, c))):
            w = 0.0; ok = True
            fs = glob.glob(d + "/*.log")
            if len(fs) < 8: continue
            for g in fs:
                s = open(g, errors="ignore").read()
                if "peak" not in s: ok = False; break
                if k == "movi": w += float(re.search(r"processing the reads: ([\d.]+) seconds", s)[1]) - float(re.search(r"bridge: setup ([\d.]+) s", s)[1])
                else: w += float(re.search(r"timing: .*\(total ([\d.]+) s\)", s)[1])
                y = re.search(r"boundaries ([\d.]+) MB, labels ([\d.]+) MB", s)
                if y: b, l = float(y[1]), float(y[2])
            if ok: v.append(w / 1.6)
        if v: rows.append((names[k], t, base[k], b, l, base[k] + Z(tags[t]), st.median(v), min(v), max(v), len(v)))
pareto = lambda r: not any(o[5] <= r[5] and o[6] <= r[6] and o != r for o in rows)
for r in rows:
    print("%-18s %-8s index %7.1f  boundaries %6.1f  labels %6.1f  whole %7.1f MB   %7.2f (%.2f-%.2f) us  n=%d %s" % (r + ("*" if pareto(r) else "",)))
print("* on the Pareto frontier (no other row is both smaller and faster)")
