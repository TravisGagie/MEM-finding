#!/usr/bin/env python3
"""Accuracy of every credit scheme on reads with k_g proportional to N_g^alpha (experiments/silva_queue.sh dial).

    dial_tables.py Q [--regions "V1_V2 V3_V4 V4_V4 V4_V5"] [--alphas "1 0 -1"] [--latex]

Q is the queue's output directory (or the collected copy, tool/results/silva_queue).  It reads
Q/dial/a<alpha>/credit/<scheme>.<region>.csv (scripts/score_silva.py, genus level) and prints:
  1. the summary: strict genus accuracy of each scheme, the mean over the regions, for each alpha;
  2. one table per alpha with each region.
Accuracy is as in scripts/silva_tables.py's tables: strict genus accuracy (a tie for the top weight counts as
wrong), the mean of the regions' accuracies; here each region's accuracy is computed from its exact count
(strict_correct / reads).  Schemes: range (lca), one, sample (cliffs), even (list), normalised (norm), proportional
(freq), the prior mixture prior<beta>, the linear mixture mix<beta> (proportional and even), and EM em<iterations>.
A missing or unusable input is printed as MISSING: <file> or BAD: <file> (<why>), never left blank; the script then
lists them and exits 1.
"""
import argparse, os, sys

NAMES = {"lca": "range (Cliffy's LCA query)", "one": "one label per MEM", "cliffs": "sample (simulated cliffs)",
         "list": "even", "norm": "normalised (c_g / N_g)", "freq": "proportional"}
FILES = ["lca", "one", "cliffs", "list", "norm", "freq", "prior", "mix", "em"]

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


def order(labels):
    def key(l):
        for i, m in enumerate(FILES):
            if l == m: return (i, 0)
            if l.startswith(m) and l[len(m):].replace(".", "", 1).isdigit(): return (i, float(l[len(m):]))
        return (99, l)
    return sorted(labels, key=key)


def name(l):
    if l in NAMES: return NAMES[l]
    if l.startswith("prior"): return "prior mixture, beta = " + l[5:]
    if l.startswith("mix"): return "proportional/even mix, beta = " + l[3:]
    if l.startswith("em"): return "EM, %s iterations" % l[2:]
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
for m in ("prior", "mix", "em"):              # a family whose file is missing everywhere still gets a row
    if not any(l.startswith(m) and l not in NAMES for l in labels): labels.append(m)


def cell(al, r, l):
    d, msgs = acc[al][r]
    if l in d: return "%.2f%%" % (100 * d[l])
    fam = next((m for m in FILES if l == m or (l.startswith(m) and m in ("prior", "mix", "em"))), l)
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
for al in alphas:
    show("alpha = %s, each region" % al, ["credit scheme", "mean"] + regions,
         [[name(l), mean(al, l)] + [cell(al, r, l) for r in regions] for l in labels])
if missing:
    print("== %d inputs missing or unreadable:" % len(missing))
    for m in missing: print("  " + m)
    sys.exit(1)
print("== every input was found")
