#!/usr/bin/env python3
"""Old and new binaries in turn on find --preset silva, one thread, pinned (tool/results/silva/README.md)."""
import glob, re, statistics as st, sys, os
W = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/research/katka/work/regress")
for v in ("old", "new"):
    m = []; t = []; pk = 0
    for d in sorted(glob.glob("%s/%s.r?" % (W, v))):
        w = n = tt = 0
        for f in glob.glob(d + "/*.log"):
            s = open(f).read(); x = re.search(r"timing: (\d+) reads, 1 threads; mean per-read time ([\d.]+) us; total time / reads ([\d.]+) us", s)
            w += int(x[1]) * float(x[2]); tt += int(x[1]) * float(x[3]); n += int(x[1])
            p = re.search(r"peak (\d+) KB", s); pk = max(pk, int(p[1]) if p else 0)
        m.append(w / n); t.append(tt / n)
    print("%-4s mean per read %.2f us (%.2f-%.2f), total/reads %.2f, peak %.2f GB, n=%d; repeats %s" % (v, st.median(m), min(m), max(m), st.median(t), pk * 1024 / 1e9, len(m), " ".join("%.2f" % x for x in m)))
