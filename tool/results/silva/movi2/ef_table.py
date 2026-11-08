#!/usr/bin/env python3
"""Summarises movi12.sh's logs (t5/<config>.r<rep>/<file>.log): µs per read over the 8 files, one thread, setup excluded."""
import glob, os, re, statistics as st, sys
D = os.path.expanduser("~/research/katka/work/pc"); M = D + "/movi"
Z = lambda f: os.path.getsize(f) / 1e6
base = {"csa": Z(D + "/bac.csa"), "ef": Z(D + "/bac.efcsa"), "movi": Z(M + "/idx/index.movi") + Z(M + "/idx/ftab.10.bin") + 58.518}
tags = {"full": D + "/bac.tag", "grammar": D + "/g.tag", "s4": D + "/bac.s4.tag"}
names = {"csa": "explicit RLCSA", "ef": "Elias-Fano RLCSA", "movi": "Movi 2, row list"}
for k in ("csa", "ef", "movi"):
    for t in ("full", "grammar", "s4"):
        c = "%s_%s" % (k, t); v = []; b = l = None
        for d in sorted(glob.glob("%s/t5/%s.r?" % (M, c))):
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
        if v:
            print("%-18s %-8s index %7.1f  boundaries %6.1f  labels %6.1f  whole %7.1f MB   %7.2f (%.2f-%.2f) us  n=%d" % (names[k], t, base[k], b, l, base[k] + Z(tags[t]), st.median(v), min(v), max(v), len(v)))
