#!/usr/bin/env python3
"""The paper's SILVA tables from the desktop's queue outputs (tool/NOTES.md, "SILVA queue: runner and tables").

    silva_tables.py Q [--regions "V1_V2 V3_V4 V4_V4 V4_V5"] [--latex]

Q is the output directory of experiments/silva_queue.sh (default tool/results/silva_queue).  The script prints,
in order:
  1. tab:credit, SILVA column: genus accuracy of each credit scheme on the same MEMs (L = 30), strict, mean over
     the regions                                                                  (step credit)
  2. tab:others, our rows: space, time per read (one thread), reads/s (all threads), genus accuracy, for the
     default configuration, grammar-compressed runs (built-in RePair and bigrepair) and minimizer digests
                                                                    (steps timing, grammar, digests)
  3. tab:digest: space and time per read for L = 15..50, strict accuracy rows      (step digests)
  4. B_lambda and the bound of Lemma "runs", rho <= (r + sigma + 2d)(lambda + 1) + B_lambda + 1     (step blam)
  5. grammar sizes: the built-in RePair against bigrepair, and whether the answers are byte-identical
                                                                                   (step grammar)
  6. the uniform-abundance simulation: reads per genus and genus accuracy per region   (step uniform)
Every number comes from a file, and every file that is missing or unreadable is named in the table, as
"MISSING: <file> (step <step>)", never left blank; a table with missing inputs ends with a list of them.
--latex also prints each table's rows in the paper's LaTeX layout (a missing value appears as \\missing{}).

Conventions (the paper's): accuracy is strict genus accuracy (a tie for the top weight counts as wrong) unless
the column says otherwise, as a mean over the regions of the per-region accuracies; time per read is the mean
per-read time with one thread over every read of every region and mate (the `timing:` line), and the total time
over the number of reads is shown beside it; reads per second is all the reads over the summed wall time of the
multithreaded runs; space is the files' total size in GB (10^9 bytes), as listed in Q/sizes.txt.
"""
import argparse, glob, math, os, re, sys

MODES = [("lca", "range (Cliffy's LCA query)"), ("one", "one label per MEM (drawn by counts)"),
         ("cliffs", "sample (simulated cliffs)"), ("list", "even"),
         ("norm", "proportional, normalised by taxon size"), ("freq", "proportional")]


class Missing(Exception):
    pass


class Ctx:
    def __init__(self, q, regions):
        self.q, self.regions, self.missing = q, regions, []

    def path(self, *p):
        return os.path.join(self.q, *p)

    def miss(self, f, step, why="missing"):
        rel = os.path.relpath(f, self.q)
        msg = "MISSING: %s (step %s)" % (rel, step) if why == "missing" else "BAD: %s (%s; step %s)" % (rel, why, step)
        if msg not in self.missing:
            self.missing.append(msg)
        return msg

    def read(self, f, step):
        if not os.path.isfile(f) or os.path.getsize(f) == 0:
            raise Missing(self.miss(f, step))
        return open(f).read()

    # accuracy at a level from a score_silva.py csv
    def score(self, f, step, level="genus"):
        txt = self.read(f, step)
        for line in txt.splitlines():
            t = line.split(",")
            if len(t) > 3 and t[0] == "katka" and t[2] == level:
                kv = dict(x.split("=", 1) for x in t[4:] if "=" in x)
                return {"acc": float(t[3]), "strict": float(kv["strict"]), "theirs": float(kv["theirs"]),
                        "reads": int(kv["reads"]), "unclassified": int(kv["unclassified"])}
        raise Missing(self.miss(f, step, "no %s line" % level))

    def mean_score(self, pattern, step, key="strict"):
        """pattern with {r}: the mean over the regions"""
        v = [self.score(self.path(pattern.format(r=r)), step)[key] for r in self.regions]
        return sum(v) / len(v)

    # the timing protocol: every region and mate, one thread and all threads
    def timing(self, d, name, step):
        one = {"reads": 0, "us": 0.0, "total": 0.0}
        many = {"reads": 0, "total": 0.0, "threads": None}
        for r in self.regions:
            for m in (1, 2):
                base = self.path(d, "%s.%s_%d" % (name, r, m))
                t = self.timing_line(base + ".t1.log", step)
                if t["threads"] != 1:
                    raise Missing(self.miss(base + ".t1.log", step, "not a one-thread run"))
                one["reads"] += t["reads"]; one["us"] += t["mean"] * t["reads"]; one["total"] += t["total"]
                logs = [f for f in glob.glob(base + ".t*.log") if not f.endswith(".t1.log")]
                if not logs:
                    raise Missing(self.miss(base + ".t<THREADS>.log", step))
                t = self.timing_line(logs[0], step)
                many["reads"] += t["reads"]; many["total"] += t["total"]; many["threads"] = t["threads"]
        return {"mean_us": one["us"] / one["reads"], "total_us": 1e6 * one["total"] / one["reads"],
                "reads_s": many["reads"] / many["total"] if many["total"] else math.inf, "threads": many["threads"],
                "reads": one["reads"]}

    def timing_line(self, f, step):
        txt = self.read(f, step)
        m = re.search(r"timing: (\d+) reads, (\d+) threads; mean per-read time ([\d.]+) us; "
                      r"total time / reads ([\d.]+) us \(total ([\d.]+) s\)", txt)
        if not m:
            raise Missing(self.miss(f, step, "no timing line"))
        if not re.search(r"exit 0\s*$", txt):
            raise Missing(self.miss(f, step, "the run did not exit 0"))
        return {"reads": int(m.group(1)), "threads": int(m.group(2)), "mean": float(m.group(3)),
                "total": float(m.group(5))}

    def size(self, name, step):
        f = self.path("sizes.txt")
        for line in self.read(f, step).splitlines():
            t = line.split("\t")
            if t[0] == name:
                return int(t[1]) / 1e9
        raise Missing(self.miss(f, step, "no line for " + name))


def cell(fn, fmt):
    """fn's value formatted, or the MISSING message"""
    try:
        return fmt(fn())
    except Missing as e:
        return str(e)


def gb(x):
    return "%.2f" % x if x >= 0.1 or x == 0 else "%.2g" % x


def pct(x):
    return "%.2f%%" % (100 * x)


def show(title, header, rows, latex=None, note=None):
    print("== " + title)
    if note:
        print(note)
    w = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    for k, r in enumerate([header] + rows):
        print("  " + "  ".join(str(x).ljust(w[i]) if i == 0 else str(x).rjust(w[i]) for i, x in enumerate(r)))
        if k == 0:
            print("  " + "  ".join("-" * x for x in w))
    if latex is not None:
        print("  LaTeX:")
        for r in latex:
            print("    " + " & ".join(("\\missing{}" if str(x).startswith(("MISSING", "BAD")) else str(x)) for x in r) + " \\\\")
    print()


def tex_pct(s):
    return s.replace("%", "\\%")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("q", nargs="?", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "silva_queue"))
    ap.add_argument("--regions", default="V1_V2 V3_V4 V4_V4 V4_V5")
    ap.add_argument("--latex", action="store_true")
    a = ap.parse_args()
    c = Ctx(os.path.abspath(a.q), a.regions.split())
    if not os.path.isdir(c.q):
        sys.exit("no directory %s: run experiments/silva_queue.sh first" % c.q)
    L = a.latex
    print("SILVA tables from %s (regions %s)\n" % (c.q, ", ".join(c.regions)))

    # 1. tab:credit
    rows = []
    for mode, label in MODES:
        rows.append([label, cell(lambda: c.mean_score("credit/%s.{r}.csv" % mode, "credit"), pct)]
                    + [cell(lambda r=r: c.score(c.path("credit/%s.%s.csv" % (mode, r)), "credit")["strict"], pct) for r in c.regions])
    show("tab:credit, SILVA (genus, strict, L = 30; mean over the regions, then each region)",
         ["credit scheme", "mean"] + c.regions, rows, [[r[0], "", tex_pct(r[1])] for r in rows] if L else None,
         "  (the Zymo column comes from the Zymo runs, not from this queue)")

    # 2. tab:others
    def others(label, tdir, tname, tstep, size, sstep, acc_pat, astep):
        t = None
        try:
            t = c.timing(tdir, tname, tstep)
        except Missing as e:
            t = str(e)
        tm = (lambda f: t if isinstance(t, str) else f(t))
        return [label, cell(lambda: c.size(size, sstep), gb),
                tm(lambda t: "%.1f" % t["mean_us"]), tm(lambda t: "%.1f" % t["total_us"]),
                tm(lambda t: "{:,}".format(int(round(t["reads_s"], -2)))), tm(lambda t: str(t["threads"])),
                cell(lambda: c.mean_score(acc_pat, astep), pct), cell(lambda: c.mean_score(acc_pat, astep, "acc"), pct)]
    rows = [others("long MEMs, proportional credit", "timing", "default", "timing", "default", "timing", "credit/freq.{r}.csv", "credit"),
            others("  grammar-compressed runs (built-in RePair)", "timing", "grammar", "grammar", "grammar", "grammar", "credit/freq.grammar.{r}.csv", "grammar"),
            others("  grammar-compressed runs (bigrepair)", "timing", "bigrepair", "grammar", "bigrepair", "grammar", "credit/freq.bigrepair.{r}.csv", "grammar"),
            others("  minimizer digests (L = 30)", "digests", "csa.L30", "digests", "digests.csa", "digests", "digests/freq.L30.{r}.csv", "digests")]
    show("tab:others, our rows (time per read: mean per-read time with one thread; total/reads beside it)",
         ["classifier", "space (GB)", "time (us)", "total/reads (us)", "reads/s", "threads", "strict acc.", "acc. (ties correct)"], rows,
         [[r[0].strip(), r[1], r[2], r[4].replace(",", "\\,"), tex_pct(r[6])] for r in rows] if L else None)

    # 3. tab:digest
    Ls = [15, 20, 30, 40, 50]
    def trow(label, name, size, Lset, step="digests"):
        out = [label, cell(lambda: c.size(size, step), gb)]
        for l in Ls:
            out.append(cell(lambda l=l: c.timing("digests", name % l, step)["mean_us"], lambda x: "%.1f" % x) if l in Lset else "")
        return out
    def arow(label, pat, Lset):
        return [label, ""] + [cell(lambda l=l: c.mean_score(pat % l, "digests"), pct) if l in Lset else "" for l in Ls]
    rows = [arow("digests: strict accuracy", "digests/freq.L%d.{r}.csv", Ls),
            trow("  explicit RLCSA, tags", "csa.L%d", "digests.csa", Ls),
            trow("  Elias-Fano RLCSA, tags", "efcsa.L%d", "digests.efcsa", Ls),
            trow("  RLBWT, tags sampled with s_r = 4", "rlbwt.s4.L%d", "digests.rlbwt.s4", [30, 40, 50]),
            arow("undigested: strict accuracy", "digests/freq.undigested.L%d.{r}.csv", [15, 20, 30, 40]),
            trow("  default configuration", "undigested.L%d", "default", [15, 20, 30, 40])]
    same = c.path("digests/identical.txt")
    note = "  (the Elias-Fano and sampled rows answer the same MEMs as the explicit row: %s)" % (
        cell(lambda: c.read(same, "digests").strip(), str))
    show("tab:digest (space in GB, then the mean per-read time in us with one thread, for L = 15 .. 50)",
         ["configuration", "space"] + ["L = %d" % l for l in Ls], rows,
         [[r[0].strip()] + [tex_pct(x) for x in r[1:]] for r in rows] if L else None, note)

    # 4. B_lambda
    f = c.path("blam/blam.log")
    try:
        txt = c.read(f, "blam")
        m = re.search(r"tag runs rho = (\d+); r = (\d+), sigma = (\d+), d = (\d+) segments", txt)
        if not m:
            raise Missing(c.miss(f, "blam", "no rho line"))
        rho, r, sig, d = map(int, m.groups())
        rows = []
        for lam, b, bound in re.findall(r"lambda = (\d+): B_lambda = (\d+), bound (\d+)", txt):
            lam, b, bound = int(lam), int(b), int(bound)
            check = (r + sig + 2 * d) * (lam + 1) + b + 1
            if check != bound:
                raise Missing(c.miss(f, "blam", "bound for lambda = %d is %d, not %d" % (lam, bound, check)))
            rows.append([str(lam), "{:,}".format(b), "{:,}".format(bound), "%.3f" % (rho / bound), "yes" if rho <= bound else "NO"])
        if not rows:
            raise Missing(c.miss(f, "blam", "no lambda lines"))
        show("B_lambda and Lemma \"runs\": rho <= (r + sigma + 2d)(lambda + 1) + B_lambda + 1",
             ["lambda", "B_lambda", "bound", "rho / bound", "holds"], rows, None,
             "  rho = {:,} tag runs; r = {:,}; sigma = {}; d = {:,} segments".format(rho, r, sig, d))
    except Missing as e:
        show("B_lambda and Lemma \"runs\"", ["lambda", "B_lambda"], [[str(e), ""]])

    # 5. grammar sizes
    f = c.path("grammar/summary.txt")
    try:
        txt = c.read(f, "grammar")
        def grab(rx, what):
            m = re.search(rx, txt)
            if not m:
                raise Missing(c.miss(f, "grammar", "no " + what))
            return m
        runs = grab(r"tag runs: (\d+) runs over (\d+) rows", "tag-run count")
        g = grab(r"1\. built-in RePair: tag array: .*?, ([\d.]+) MB; wall ([\d.]+) s, peak ([\d.]+) GB", "built-in RePair line")
        dnf = re.search(r"2\. bigrepair -i: DID NOT FINISH.*", txt)
        if dnf:
            c.miss(f, "grammar", "bigrepair did not finish: " + dnf.group(0)[len("2. bigrepair -i: "):])
        bg = None if dnf else grab(r"2\. bigrepair -i: \.R (\d+) bytes, \.C (\d+) bytes; wall ([\d.]+) s, peak ([\d.]+) GB", "bigrepair line")
        bt = None if dnf else grab(r"3\. from bigrepair's grammar: tag array: .*?, ([\d.]+) MB; wall ([\d.]+) s, peak ([\d.]+) GB", "bigrepair tag array line")
        fl = grab(r"tag array files: plain (\d+) bytes \(\+ RMQ (\d+)\), built-in RePair (\d+), bigrepair (\d+|none)", "file sizes")
        ident = re.findall(r"find -(\w): outputs byte-identical with the (?:three )?tag arrays[^:]*: (\w+)", txt)
        mb = lambda x: gb(int(x) / 1e6)
        rows = [["plain (+ RMQ)", "%s (+ %s)" % (mb(fl.group(1)), mb(fl.group(2))), "", ""],
                ["built-in RePair", mb(fl.group(3)), "%s s" % g.group(2), "%s GB" % g.group(3)],
                ["bigrepair", mb(fl.group(4)), "%.1f s" % (float(bg.group(3)) + float(bt.group(2))),
                 "%s GB" % max(bg.group(4), bt.group(3), key=float)] if bg else
                ["bigrepair", "BAD: grammar/summary.txt (bigrepair did not finish; step grammar)", "", ""]]
        show("Grammar-compressed tags (file sizes in MB; build time and peak memory)",
             ["tag array", "file (MB)", "build", "peak"], rows, None,
             "  {:,} tag runs over {:,} rows; bigrepair's grammar: .R {} MB, .C {} MB; answers byte-identical: {}".format(
                 int(runs.group(1)), int(runs.group(2)), mb(bg.group(1)) if bg else "-", mb(bg.group(2)) if bg else "-",
                 ", ".join("-%s %s" % x for x in ident) if len(ident) == 2 else c.miss(f, "grammar", "the find comparison did not finish")))
    except Missing as e:
        show("Grammar-compressed tags", ["tag array", "file (MB)"], [[str(e), ""]])

    # 6. the uniform simulation
    rows = []
    for r in c.regions:
        f = c.path("uniform/%s.counts.txt" % r)
        def counts(f=f):
            v = [int(x) for x in c.read(f, "uniform").split()]
            nz = [x for x in v if x > 0]
            if not nz:
                raise Missing(c.miss(f, "uniform", "no reads"))
            dr = c.path("uniform/%s.dropped.txt" % r)
            nd = sum(1 for l in open(dr) if l.strip()) if os.path.exists(dr) else None
            return "{:,} pairs, {:,} of {:,} genera, {}..{} each; {}".format(sum(v), len(nz), len(v), min(nz), max(nz),
                "%d dropped (too few reads)" % nd if nd is not None else c.miss(dr, "uniform"))
        sc = lambda k, r=r: c.score(c.path("uniform/%s.score.csv" % r), "uniform")[k]
        rows.append([r, cell(counts, str), cell(lambda: sc("strict"), pct), cell(lambda: sc("acc"), pct)])
    rows.append(["mean", "", cell(lambda: c.mean_score("uniform/{r}.score.csv", "uniform"), pct),
                 cell(lambda: c.mean_score("uniform/{r}.score.csv", "uniform", "acc"), pct)])
    show("Uniform-abundance simulation (proportional credit, default configuration, genus)",
         ["region", "reads", "strict acc.", "acc. (ties correct)"], rows)

    if c.missing:
        print("== %d inputs missing or unreadable:" % len(c.missing))
        for m in c.missing:
            print("  " + m)
        sys.exit(1)
    print("== every input was found")


if __name__ == "__main__":
    main()
