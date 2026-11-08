#!/usr/bin/env python3
"""The time-space table from timing14.sh's logs (T14/<config>.p<threads>.r<rep>/<file>.{all,one}.time, GNU time's
"wall_user_sys_peakKB"): every tool timed the same way, wall time on all the reads minus wall time on a single read
(loading and setup), summed over the files.  One thread (pinned): microseconds per read, median of the repeats (range).
12 threads: reads per second, median.  Sizes in MB (10^6 bytes): the index in memory with its lookup table (RLCSAs: the 17.8 MB
10-mer table built at load; Movi 2: its ftab and the 58.5 MB row list), plus the tag array; Kraken 2: its .k2d files;
Tagger: its index files.  Also: KATKA's own timing (find's "timing:" total, setup excluded) for the RLCSA rows."""
import glob, os, re, statistics as st, sys
H = os.path.expanduser("~/research/katka/work"); P = H + "/pc"; M = P + "/movi"; T14 = sys.argv[1] if len(sys.argv) > 1 else H + "/t14"
Z = lambda *f: sum(os.path.getsize(x) for x in f) / 1e6
NREADS = 1_600_000
base = {"csa": Z(P + "/bac.csa") + 17.826, "ef": Z(P + "/bac.efcsa") + 8.54 + 17.826,   # + the select samples built at load (KS 4.27, FP 4.27 MB)
         "movi": Z(M + "/idx/index.movi", M + "/idx/ftab.10.bin") + 58.518}
tags = {"full": P + "/bac.tag", "grammar": P + "/g.tag", "s2": P + "/bac.s2.tag", "s4": P + "/bac.s4.tag", "s8": P + "/bac.s8.tag"}
rows = [("movi", "full"), ("movi", "grammar"), ("csa", "full"), ("csa", "s2"), ("csa", "grammar"), ("csa", "s4"), ("csa", "s8"),
        ("ef", "grammar"), ("ef", "s4"), ("ef", "s8"), ("kraken2", ""), ("tagger_L25", ""), ("tagger_L30", "")]
names = {"movi": "Movi 2, row list", "csa": "explicit RLCSA", "ef": "Elias-Fano RLCSA", "kraken2": "Kraken 2 (SILVA build)", "tagger_L25": "Tagger, L = 25", "tagger_L30": "Tagger, L = 30"}
tagname = {"full": "RLE", "grammar": "SLP", "s2": "s = 2", "s4": "s = 4", "s8": "s = 8", "": ""}
def wall(f): return float(open(f).read().split()[-1].split("_")[0])
def run(cfg, p):
    out = []; internal = []
    for d in sorted(glob.glob("%s/%s.p%d.r?" % (T14, cfg, p))):
        A = sorted(glob.glob(d + "/*.all.time"))
        need = 4 if cfg == "kraken2" else 8
        if len(A) < need or len(glob.glob(d + "/*.one.time")) < need: continue
        out.append(sum(wall(a) - wall(a.replace(".all.", ".one.")) for a in A))
        if cfg.startswith(("csa", "ef")) and p == 1:
            internal.append(sum(float(re.search(r"\(total ([\d.]+) s\)", open(a.replace(".time", ".log")).read())[1]) for a in A))
    return out, internal
print("%-24s %-6s %9s  %-26s %-16s %s" % ("index", "tags", "MB", "us/read, 1 thread", "reads/s, 12 thr.", "KATKA's own timing, 1 thread"))
for k, t in rows:
    cfg = "%s_%s" % (k, t) if t else k
    if k in base: size = base[k] + Z(tags[t])
    elif k == "kraken2": size = Z(*glob.glob(H + "/k2/std/*.k2d"))
    else: size = Z(*[f for f in glob.glob(H + "/tagger/full/silva*") if os.path.isfile(f)])
    one, internal = run(cfg, 1); twelve, _ = run(cfg, 12)
    s1 = "%7.2f (%.2f-%.2f) n=%d" % (st.median(one) / NREADS * 1e6, min(one) / NREADS * 1e6, max(one) / NREADS * 1e6, len(one)) if one else "-"
    s12 = "%8.0f n=%d" % (NREADS / st.median(twelve), len(twelve)) if twelve else "-"
    si = "%7.2f us" % (st.median(internal) / NREADS * 1e6) if internal else ""
    print("%-24s %-6s %9.1f  %-26s %-16s %s" % (names[k], tagname[t], size, s1, s12, si))
