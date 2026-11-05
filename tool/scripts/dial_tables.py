#!/usr/bin/env python3
"""Accuracy of every credit scheme on reads with k_g proportional to N_g^alpha (experiments/silva_queue.sh dial).

    dial_tables.py Q [--regions "V1_V2 V3_V4 V4_V4 V4_V5"] [--alphas "1 0 -1"] [--latex]

Q is the queue's output directory (or the collected copy, tool/results/silva_queue).  It reads
Q/dial/a<alpha>/credit/<scheme>.<region>.csv (scripts/score_silva.py, genus level) and prints:
  1. the summary: strict genus accuracy of each scheme, the mean over the regions, for each alpha;
  2. the held-out choice of epsilon for naive Bayes with the flat prior: the epsilon with the best mean at alpha = 0,
     its accuracy at each alpha, and the best epsilon for each alpha for comparison;
  3. one table per alpha with each region.
Accuracy is as in scripts/silva_tables.py's tables: strict genus accuracy (a tie for the top weight counts as
wrong), the mean of the regions' accuracies; here each region's accuracy is computed from its exact count
(strict_correct / reads).  Schemes: range (lca), one, sample (cliffs), even (list), normalised (norm), proportional
(freq), the prior mixture prior<beta>, the linear mixture mix<beta> (proportional and even), EM em<iterations>, and
naive Bayes over the MEMs: with the sample's own prior (bayes, the reference), a flat prior (bayesflat, with
bayesflat<epsilon> for the sweep of the missing-match likelihood) and the database prior (bayesdb).
A missing or unusable input is printed as MISSING: <file> or BAD: <file> (<why>), never left blank; the script then
lists them and exits 1.
"""
import argparse, os, sys

NAMES = {"lca": "range (Cliffy's LCA query)", "one": "one label per MEM", "cliffs": "sample (simulated cliffs)",
         "list": "even", "norm": "normalised (c_g / N_g)", "freq": "proportional",
         "bayes": "naive Bayes, the sample's prior (reference; epsilon = 1e-06)", "bayesflat": "naive Bayes, flat prior, epsilon = 1e-06",
         "bayesdb": "naive Bayes, database prior, epsilon = 1e-06",
         "bayesauto": "naive Bayes, mixture with epsilon fitted to the reads"}
FILES = ["lca", "one", "cliffs", "list", "norm", "freq", "prior", "mix", "em", "bayes", "bayesflat", "bayesdb", "bayesauto"]
FAMILIES = ("prior", "mix", "em", "bayesflat")       # files with several lines (one per beta, iteration count or epsilon)

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("q", nargs="?", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "silva_queue"))
ap.add_argument("--regions", default="V1_V2 V3_V4 V4_V4 V4_V5")
ap.add_argument("--alphas", default="1 0 -1")
ap.add_argument("--latex", action="store_true")
a = ap.parse_args()
Q = os.path.abspath(a.q); regions = a.regions.split(); alphas = a.alphas.split()
missing = []


def miss(f, why=None):
    m = "MISSING: %s (step dial)" % os.path.relpath(f, Q) if why is None else "BAD: %s (%s; step dial)" % (os.path.relpath(f, Q), why)
    if m not in missing: missing.append(m)
    return m


def load(al, mode, r):
    """{label: strict accuracy} at genus level from one csv, or a message"""
    f = os.path.join(Q, "dial", "a" + al, "credit", "%s.%s.csv" % (mode, r))
    if not os.path.isfile(f) or os.path.getsize(f) == 0:
        return miss(f)
    out = {}
    for line in open(f):
        t = line.strip().split(",")
        if len(t) < 4 or t[0] != "katka" or t[2] != "genus": continue
        kv = dict(x.split("=", 1) for x in t[4:] if "=" in x)
        if "strict_correct" not in kv:
            return miss(f, "no strict_correct field: scored by an older score_silva.py")
        out[t[1]] = int(kv["strict_correct"]) / int(kv["reads"])
    return out or miss(f, "no genus lines")


def family(l):
    """the file a label comes from"""
    if l in FILES: return l
    for m in sorted(FAMILIES, key=len, reverse=True):
        if l.startswith(m):
            try: float(l[len(m):]); return m
            except ValueError: pass
    return None


def order(labels):
    def key(l):
        m = family(l)
        if m is None: return (99, 0, l)
        v = l[len(m):]
        # epsilon sweeps from the largest epsilon down; the default (no suffix, 1e-6) in its place
        val = float(v) if v else (1e-6 if m == "bayesflat" else 0)
        return (FILES.index(m), -val if m == "bayesflat" else val, l)
    return sorted(labels, key=key)


def name(l):
    if l in NAMES: return NAMES[l]
    if l.startswith("prior"): return "prior mixture, beta = " + l[5:]
    if l.startswith("mix"): return "proportional/even mix, beta = " + l[3:]
    if l.startswith("em"): return "EM, %s iterations" % l[2:]
    if l.startswith("bayesflat"): return "naive Bayes, flat prior, epsilon = %s" % (l[9:] or "1e-06")
    return l


# acc[al][r]: {label: accuracy} or a message per file
acc = {}
for al in alphas:
    acc[al] = {}
    for r in regions:
        d = {}; msgs = {}
        for m in FILES:
            x = load(al, m, r)
            if isinstance(x, str): msgs[m] = x
            else: d.update(x)
        acc[al][r] = (d, msgs)
labels = order({l for al in alphas for r in regions for l in acc[al][r][0]} | set(FILES[:6]))
for m in FAMILIES:                            # a family whose file is missing everywhere still gets a row
    if not any(family(l) == m and l != m for l in labels) and m not in labels: labels.append(m)


def cell(al, r, l):
    d, msgs = acc[al][r]
    if l in d: return "%.2f%%" % (100 * d[l])
    fam = family(l) or l
    return msgs.get(fam) or miss(os.path.join(Q, "dial", "a" + al, "credit", "%s.%s.csv" % (fam, r)), "no line " + l)


def mean(al, l):
    v = [acc[al][r][0].get(l) for r in regions]
    if any(x is None for x in v): return cell(al, regions[[x is None for x in v].index(True)], l)
    return "%.2f%%" % (100 * sum(v) / len(v))


def show(title, header, rows):
    print("== " + title)
    w = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    for k, r in enumerate([header] + rows):
        print("  " + "  ".join(str(x).ljust(w[i]) if i == 0 else str(x).rjust(w[i]) for i, x in enumerate(r)))
        if k == 0: print("  " + "  ".join("-" * x for x in w))
    if a.latex:
        print("  LaTeX:")
        for r in rows:
            print("    " + " & ".join("\\missing{}" if str(x).startswith(("MISSING", "BAD")) else str(x).replace("%", "\\%") for x in r) + " \\\\")
    print()


print("Credit schemes on reads with k_g ~ N_g^alpha, from %s (regions %s)\n" % (Q, ", ".join(regions)))
for al in alphas:                                 # the reads behind each alpha
    info = []
    for r in regions:
        f = os.path.join(Q, "dial", "a" + al, "%s.seqtax.txt" % r)
        if os.path.isfile(f):
            c = [int(l.rstrip("\n").split("\t")[-1]) for l in open(f) if l.strip()]
            nz = [x for x in c if x]
            info.append("%s %d pairs over %d genera (%d..%d)" % (r, sum(c), len(nz), min(nz), max(nz)) if nz else r + " no reads")
        else:
            info.append(miss(f))
    print("alpha = %s: %s" % (al, "; ".join(info)))
print()
show("Strict genus accuracy, mean over the regions", ["credit scheme"] + ["alpha = %s" % al for al in alphas],
     [[name(l)] + [mean(al, l) for al in alphas] for l in labels])
# the held-out choice of epsilon for naive Bayes with the flat prior: chosen on alpha = 0, reported on the others
eps = [l for l in labels if family(l) == "bayesflat"]
def mval(al, l):
    v = [acc[al][r][0].get(l) for r in regions]
    return None if any(x is None for x in v) else sum(v) / len(v)
if "0" in alphas and eps:
    m0 = {l: mval("0", l) for l in eps}
    if any(v is None for v in m0.values()):
        print("== epsilon chosen on alpha = 0: %s\n" % mean("0", next(l for l in eps if m0[l] is None)))
    else:
        best = max(eps, key=lambda l: (m0[l], -eps.index(l)))      # ties: the larger epsilon (listed first)
        rows = []
        for al in alphas:
            vals = {l: mval(al, l) for l in eps}
            own = max((l for l in eps if vals[l] is not None), key=lambda l: vals[l], default=None)
            rows.append(["alpha = %s" % al, mean(al, best), name(own).split("epsilon = ")[-1] if own else "-",
                         mean(al, own) if own else "-"])
        show("Naive Bayes, flat prior: epsilon chosen on alpha = 0 (%s), and the best epsilon for each alpha"
             % name(best).split("epsilon = ")[-1], ["", "chosen epsilon", "best epsilon", "its accuracy"], rows)
# the epsilon bayesauto fitted, per alpha and region (from its log)
rows = []
for al in alphas:
    row = ["alpha = %s" % al]
    for r in regions:
        f = os.path.join(Q, "dial", "a" + al, "credit", "bayesauto.%s.csv.log" % r)
        m = None
        if os.path.isfile(f):
            for line in open(f):
                if line.startswith("bayesauto: epsilon = "): m = line.split("=")[1].split("(")[0].strip() + (" (edge)" if "edge of the grid" in line else "")
        row.append(m if m else miss(f))
    rows.append(row)
show("Naive Bayes (mixture): the epsilon fitted to each read set", [""] + regions, rows)
for al in alphas:
    show("alpha = %s, each region" % al, ["credit scheme", "mean"] + regions,
         [[name(l), mean(al, l)] + [cell(al, r, l) for r in regions] for l in labels])
if missing:
    print("== %d inputs missing or unreadable:" % len(missing))
    for m in missing: print("  " + m)
    sys.exit(1)
print("== every input was found")
