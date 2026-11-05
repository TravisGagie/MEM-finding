#!/usr/bin/env python3
# Scores SILVA listings (katka-classify output for each mate) at every rank, under the credit schemes of Table 1:
#   --mode freq = proportional, list = even, norm = normalized, one = one draw, cliffs = sample, lca = range,
#   longest = only the longest MEMs of each mate, with proportional credit.  Needs cliffy-experiments (--exp-src).
# (From silva/rz_score_schemes.py in the original repository; only the names have changed.)
# Fast, equivalent re-implementation of "method0" (Cliffy's LCA query) from cliffy-experiments'
# classify_silva_readset.py.  Their version deep-copies a weight table with one entry per clade for
# every read and sends it to worker processes with each read, which takes hours for 1M read pairs on a
# desktop.  This one loads the data with their own classes and applies the same rule with a sparse
# weight table:
#   each match of length L in [left doc, right doc] adds L/(right-left+1) to every document's clade
#   in that range; the read is correct if its true clade is among the clades of maximum weight
#   (so a read with no weight at all counts as correct unless its clade is "uncultured", as in their code); reads whose true clade is
#   missing or ambiguous at a level count as correct, as in their code.
# Schemes added for samples that don't follow the database's composition (each is computed from the same {doc:count}
# listings; N_g = sequences of genus g (--nseq), N = sum of N_g, G = number of genera; every value of beta or of
# the number of iterations is a line of its own, labelled e.g. prior0.25, mix0.5, em20):
#   --mode prior: credit to g proportional to (c_g / N_g) * ((1 - beta) N_g / N + beta / G), for each --beta;
#                 beta = 0 is proportional credit (bit for bit), beta = 1 normalised credit (c_g / N_g)
#   --mode mix:   (1 - beta) * proportional + beta * even, per match, for each --beta (beta = 0 and 1 bit for bit)
#   --mode em:    credit proportional to (c_g / N_g) * pi_g, where pi starts as the database prior N_g / N and is
#                 re-estimated as each genus's share of the total credit over all reads and both mates, for each
#                 number of iterations in --em-iters (0 iterations is proportional credit, up to rounding)
# Output: katka,<mode>,<level>,<accuracy>,correct=..,unclassified=..,wrong=..,reads=..,ambiguous_truth=..,theirs=..,strict=..,strict_correct=..
# (accuracy counts unclassified reads as wrong; theirs = accuracy under their rule for unclassified reads;
#  strict = accuracy when a tie for the maximum weight counts as wrong; all three count ties as in their code otherwise).
# --mode lca: a feature {l,r} spreads its length over documents l..r (method0); --mode list: a feature
# {d1..dk} gives length/k to each listed document (the method2 weighting, |M|/|L|); in list mode a feature
# <l,r> (hybrid: answered by its LCA) is weighted as in lca mode.
import sys, os, argparse, importlib.util, random
import multiprocessing
# fork, so workers share the loaded read set (Python 3.14 defaults to forkserver on Linux)
Pool = multiprocessing.get_context("fork").Pool

p = argparse.ArgumentParser()
p.add_argument("--mate1-listings", required=True); p.add_argument("--mate2-listings", required=True)
p.add_argument("--doc-id-to-traversal", required=True); p.add_argument("--silva-tax-ranks", required=True)
p.add_argument("--readset-truthset", required=True); p.add_argument("--trav-to-length", required=True)
p.add_argument("--output", required=True); p.add_argument("--exp-src", required=True)
p.add_argument("--threads", type=int, default=os.cpu_count())
p.add_argument("--mode", choices=["lca", "list", "freq", "norm", "one", "cliffs", "longest", "prior", "mix", "em"], default="lca")
p.add_argument("--beta", default="0,0.1,0.25,0.5,0.75,1", help="--mode prior or mix: the values of beta (comma-separated)")
p.add_argument("--em-iters", default="5,20", help="--mode em: the numbers of EM iterations to report (comma-separated)")
p.add_argument("--seed", type=int, default=1)
p.add_argument("--levels", default="genus,family,order,class,phylum,domain", help="the levels to score (comma-separated)")
p.add_argument("--nseq", default=None, help="sequences per document (one per line), for --mode norm, prior and em")
a = p.parse_args()
spec = importlib.util.spec_from_file_location("csr", os.path.join(a.exp_src, "classify_silva_readset.py"))
csr = importlib.util.module_from_spec(spec); spec.loader.exec_module(csr)
rs = csr.DocProfReadSet(a.mate1_listings, a.mate2_listings, a.doc_id_to_traversal, a.silva_tax_ranks,
                        a.readset_truthset, a.trav_to_length)

def parse(listing):
    t = listing.split(); out = []
    for i in range(0, len(t), 2):
        s, e = t[i][1:-1].split(",")
        d, c = [], []
        for x in t[i + 1][1:-1].split(","):
            if ":" in x: a, b = x.split(":"); d.append(int(a)); c.append(int(b))
            else: d.append(int(x)); c.append(1)
        out.append((int(e) - int(s) + 1, d, t[i + 1][0] == "<", c))   # <l,r>: an LCA answer in a hybrid listing
    return out
NSEQ = [int(l) for l in open(a.nseq)] if a.nseq else None
if a.mode in ("norm", "prior", "em") and NSEQ is None:
    sys.exit("--mode %s needs --nseq" % a.mode)
NTOT, NGEN = (sum(NSEQ), len(NSEQ)) if NSEQ else (0, 0)
PI0 = [n / NTOT for n in NSEQ] if NSEQ else None      # the database prior, EM's start

G = {}
def work(rng):
    lo, hi = rng; keys = G["keys"]; d2c = G["d2c"]; labels = G["labels"]; level = G["level"]
    tp = uncl = uncl_theirs = amb = strict = 0
    listmode = G["mode"] in ("list", "freq", "norm", "one", "cliffs", "longest", "prior", "mix", "em")
    beta = G.get("beta"); pi = G.get("pi")
    rnd = random.Random(a.seed * 1000003 + lo)
    for r in rs.read_list[lo:hi]:
        correct = r.correct_taxa.get_certain_level(level)
        n = labels.get(correct, 0)
        if n != 1: tp += 1; amb += 1; strict += 1; continue      # truth missing/ambiguous at this level: correct, as in their code
        acc = {}
        for listing in (r.mate1_listings, r.mate2_listings):
            feats = parse(listing)
            if G["mode"] == "longest" and feats:      # only the longest MEMs of each mate vote (e.g. katka-classify -L 1 -X)
                top = max(f[0] for f in feats); feats = [f for f in feats if f[0] == top]
            for L, d, rng, cnt in feats:
                if not rng and G["mode"] in ("prior", "em"):
                    # prior: (c_g / N_g) ((1 - beta) N_g / N + beta / G), written as c_g ((1 - beta) G N_g + beta N) / (G N_g)
                    # so that beta = 0 gives c_g exactly; em: (c_g / N_g) pi_g
                    # (the same weights up to a constant factor, which each match's normalisation removes; the forms
                    # are chosen so that beta = 0 is freq and beta = 1 is norm bit for bit, and so is EM's start)
                    if G["mode"] == "prior" and beta < 1:
                        f = [c * ((1 - beta) * (NGEN * NSEQ[doc]) + beta * NTOT) / (NGEN * NSEQ[doc]) for doc, c in zip(d, cnt)]
                    elif G["mode"] == "prior":
                        f = [c / (2 * NSEQ[doc]) for doc, c in zip(d, cnt)]
                    else:
                        f = [c * (pi[doc] / PI0[doc]) for doc, c in zip(d, cnt)]
                    tot = sum(f)
                    if tot <= 0: continue
                    for doc, x in zip(d, f):
                        t = d2c.get(doc)
                        if t and t in keys: acc[t] = acc.get(t, 0) + L * x / tot
                    continue
                if not rng and G["mode"] == "mix":
                    # (1 - beta) * proportional + beta * even, per match
                    tot = sum(cnt); w_even = L / len(d)
                    for doc, c in zip(d, cnt):
                        t = d2c.get(doc)
                        if t and t in keys: acc[t] = acc.get(t, 0) + ((1 - beta) * (L * c / tot) + beta * w_even)
                    continue
                if listmode and not rng and G["mode"] in ("freq", "norm", "longest"):
                    # credit proportional to the occurrences in each genus (freq), or to the fraction of the
                    # genus's sequences (both strands) containing the match (norm)
                    f = [c if G["mode"] != "norm" else c / (2 * NSEQ[doc]) for doc, c in zip(d, cnt)]
                    tot = sum(f)
                    for doc, x in zip(d, f):
                        t = d2c.get(doc)
                        if t and t in keys: acc[t] = acc.get(t, 0) + L * x / tot
                    continue
                if G["mode"] == "one" and not rng:
                    # all the credit to one genus, drawn with probability proportional to its occurrences
                    x = rnd.random() * sum(cnt); doc = d[-1]
                    for dd, c in zip(d, cnt):
                        x -= c
                        if x < 0: doc = dd; break
                    t = d2c.get(doc)
                    if t and t in keys: acc[t] = acc.get(t, 0) + L
                    continue
                if G["mode"] == "cliffs" and not rng:
                    # an even split over a cliffs-shaped sample: random values, left-to-right and right-to-left maxima
                    v = [rnd.random() for _ in d]; keep = set(); m = -1.0
                    for i, y in enumerate(v):
                        if y > m: keep.add(i); m = y
                    m = -1.0
                    for i in range(len(v) - 1, -1, -1):
                        if v[i] > m: keep.add(i); m = v[i]
                    docs = [d[i] for i in sorted(keep)]
                else:
                    docs = d if listmode and not rng else range(d[0], d[-1] + 1)
                w = L / len(docs)
                for doc in docs:
                    t = d2c.get(doc)
                    if t and t in keys: acc[t] = acc.get(t, 0) + w
        if not acc:                                  # unclassified (their code: correct unless "uncultured")
            uncl += 1; uncl_theirs += correct in keys; continue
        m = max(acc.values())
        if m <= 0 or acc.get(correct) == m:
            tp += 1
            if sum(1 for v in acc.values() if v == m) == 1: strict += 1   # the true clade alone has the maximum
    return (tp, uncl, uncl_theirs, amb, strict)

def em_pass(rng):
    """the total credit per document over reads lo..hi, under the current prior G["pi"]"""
    lo, hi = rng; pi = G["pi"]; T = [0.0] * NGEN
    for r in rs.read_list[lo:hi]:
        for listing in (r.mate1_listings, r.mate2_listings):
            for L, d, rng_, cnt in parse(listing):
                if rng_: continue
                f = [c * (pi[doc] / PI0[doc]) for doc, c in zip(d, cnt)]
                tot = sum(f)
                if tot <= 0: continue
                for doc, x in zip(d, f): T[doc] += L * x / tot
    return T

N = len(rs.read_list)
step = (N + 4 * a.threads - 1) // (4 * a.threads)
chunks = [(i, min(N, i + step)) for i in range(0, N, step)]
# the variants of the scheme, as (label, settings): one per beta, or per number of EM iterations
if a.mode in ("prior", "mix"):
    variants = [("%s%s" % (a.mode, b), {"beta": float(b)}) for b in a.beta.split(",")]
elif a.mode == "em":
    iters = sorted({int(k) for k in a.em_iters.split(",")})
    pi = list(PI0); variants = []
    for k in range(iters[-1] + 1):
        if k in iters: variants.append(("em%d" % k, {"pi": pi}))
        if k == iters[-1]: break
        G["pi"] = pi
        with Pool(a.threads) as pool:
            parts = pool.map(em_pass, chunks)
        T = [sum(x[i] for x in parts) for i in range(NGEN)]; tt = sum(T)
        pi = [x / tt for x in T]
        print("EM iteration %d: %d genera with credit" % (k + 1, sum(1 for x in pi if x > 0)), file=sys.stderr, flush=True)
else:
    variants = [(a.mode, {})]
with open(a.output, "w") as out:
  for label, settings in variants:
    G.update(settings)
    for level in a.levels.split(","):
        labels = {}
        for x in rs.silva_nodes_at_each_level[level]:
            c = x.get_certain_level(level); labels[c] = labels.get(c, 0) + 1
        G.update(mode=a.mode, level=level, labels=labels, keys={c for c in labels if c != "uncultured"},
                 d2c={i: x.get_certain_level(level) for i, x in rs.doc_id_to_full_traversal.items()})
        with Pool(a.threads) as pool:           # forked workers see G and rs without pickling them
            res = pool.map(work, chunks)
        tp, uncl, uth, amb, st = (sum(x[i] for x in res) for i in range(5))
        # accuracy counts unclassified reads as wrong; "theirs" applies their rule to unclassified reads
        line = (f"katka,{label},{level},{round(tp / N, 4)},correct={tp},unclassified={uncl},wrong={N - tp - uncl},"
                f"reads={N},ambiguous_truth={amb},theirs={round((tp + uth) / N, 4)},strict={round(st / N, 4)},strict_correct={st}")
        print(line, flush=True); out.write(line + "\n")
