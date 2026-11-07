#!/usr/bin/env python3
"""Summarises movi9.sh's logs, t3/<config>.r<rep>/<file>.log.  Time per read: the whole processing time over the 200,000 reads of each
of the 8 files, without loading; for Movi 2, its processing time minus the prototype's setup, which falls inside that clock."""
import glob, os, re, statistics as st, sys
D = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/tool/silva/pc/movi")
for c in (sys.argv[2] if len(sys.argv) > 2 else "csafull oldfull selfull csas2 olds2 sels2 csas4 olds4 sels4 csas8 olds8 sels8 olds16 csag oldg selg oldbg").split():
    t = []; pk = 0; P = S = 0; lf = 0
    for d in sorted(glob.glob("%s/t3/%s.r?" % (D, c))):
        w = 0.0
        for g in glob.glob(d + "/*.log"):
            s = open(g, errors="ignore").read()
            if "peak" not in s: w = None; break
            if "csa" in c: w += float(re.search(r"timing: .*\(total ([\d.]+) s\)", s)[1])
            else:
                w += float(re.search(r"processing the reads: ([\d.]+) seconds", s)[1]) - float(re.search(r"bridge: setup ([\d.]+) s", s)[1])
                y = re.search(r"BRIDGE points (\d+) scanned (\d+) lf_steps (\d+)", s)
                if y: P += int(y[1]); S += int(y[2]); lf += int(y[3])
            p = re.search(r"peak (\d+) KB", s); pk = max(pk, int(p[1]) if p else 0)
        if w is not None and len(glob.glob(d + "/*.log")) == 8: t.append(w / 1.6)
    if t: print("%-9s %8.2f (%.2f-%.2f) peak %.2f n=%d points %.1fM scanned/point %.2f lf_steps %.1fM" % (c, st.median(t), min(t), max(t), pk * 1024 / 1e9, len(t), P / len(t) / 1e6, S / max(P, 1), lf / len(t) / 1e6))
