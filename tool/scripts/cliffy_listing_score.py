#!/usr/bin/env python3
"""Scoring for Cliffy's approximate document listing: a fast equivalent of "method2" of cliffy-experiments'
classify_silva_readset.py (Ahmed, Boucher and Langmead's scoring), as katka/silva/rz_score2.py is of "method0"
(the LCA query).

    cliffy_listing_score.py --mate1-listings M1 --mate2-listings M2 --doc-id-to-traversal F --silva-tax-ranks F
                            --readset-truthset F --trav-to-length F --exp-src DIR --output OUT.csv [--threads T]

The listings are `cliffy run --taxcomp` output: per read, `[start,end] {d1,d2,...}` per match.  The rule, as in
their classify_method2: every clade at the level (but "uncultured") starts at weight 0; each match of length
end - start + 1 adds its full length to the clade of every listed document (once per listed document); the read
is correct if its true clade is among the clades of maximum weight (ties count as correct, and so does a read
with no match at all, unless its clade is "uncultured"); a read whose true clade is missing or ambiguous at the
level counts as correct.  The data are loaded with their own classes (--exp-src: cliffy-experiments/src), so the
files are read exactly as their script reads them.

Output, per level (genus to domain), one line:
    cliffy_listing,<level>,<accuracy>,correct=..,unclassified=..,wrong=..,reads=..,ambiguous_truth=..,theirs=..,strict=..
accuracy counts reads without a match as wrong; theirs applies their rule to them (the number their script
reports as "cliffy,method2,<level>,..."); strict also counts ties for the maximum as wrong.
"""
import argparse, importlib.util, multiprocessing, os

Pool = multiprocessing.get_context("fork").Pool

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("mate1-listings", "mate2-listings", "doc-id-to-traversal", "silva-tax-ranks", "readset-truthset", "trav-to-length", "exp-src", "output"):
    p.add_argument("--" + x, required=True)
p.add_argument("--threads", type=int, default=os.cpu_count())
a = p.parse_args()
spec = importlib.util.spec_from_file_location("csr", os.path.join(a.exp_src, "classify_silva_readset.py"))
csr = importlib.util.module_from_spec(spec); spec.loader.exec_module(csr)
rs = csr.DocProfReadSet(a.mate1_listings, a.mate2_listings, a.doc_id_to_traversal, a.silva_tax_ranks,
                        a.readset_truthset, a.trav_to_length)


def parse(listing):
    t = listing.split(); out = []
    assert len(t) % 2 == 0
    for i in range(0, len(t), 2):
        s, e = t[i][1:-1].split(",")
        out.append((int(e) - int(s) + 1, [int(x) for x in t[i + 1][1:-1].split(",")]))
    return out


G = {}


def work(rng):
    lo, hi = rng; keys = G["keys"]; d2c = G["d2c"]; labels = G["labels"]; level = G["level"]
    tp = uncl = uncl_theirs = amb = strict = 0
    for r in rs.read_list[lo:hi]:
        correct = r.correct_taxa.get_certain_level(level)
        if labels.get(correct, 0) != 1:              # truth missing or ambiguous at this level: correct, as theirs
            tp += 1; amb += 1; strict += 1; continue
        acc = {}
        for listing in (r.mate1_listings, r.mate2_listings):
            for L, docs in parse(listing):
                for doc in docs:
                    t = d2c.get(doc)
                    if t and t in keys: acc[t] = acc.get(t, 0) + L
        if not acc:                                  # no match: every clade ties at 0 (correct unless "uncultured")
            uncl += 1; uncl_theirs += correct in keys; continue
        m = max(acc.values())
        if acc.get(correct) == m:
            tp += 1
            if sum(1 for v in acc.values() if v == m) == 1: strict += 1
    return tp, uncl, uncl_theirs, amb, strict


N = len(rs.read_list)
with open(a.output, "w") as out:
    for level in ["genus", "family", "order", "class", "phylum", "domain"]:
        labels = {}
        for x in rs.silva_nodes_at_each_level[level]:
            c = x.get_certain_level(level); labels[c] = labels.get(c, 0) + 1
        G.update(level=level, labels=labels, keys={c for c in labels if c != "uncultured"},
                 d2c={i: x.get_certain_level(level) for i, x in rs.doc_id_to_full_traversal.items()})
        step = max(1, (N + 4 * a.threads - 1) // (4 * a.threads))
        with Pool(a.threads) as pool:
            res = pool.map(work, [(i, min(N, i + step)) for i in range(0, N, step)])
        tp, uncl, uth, amb, st = (sum(x[i] for x in res) for i in range(5))
        line = (f"cliffy_listing,{level},{round(tp / N, 4)},correct={tp},unclassified={uncl},wrong={N - tp - uncl},"
                f"reads={N},ambiguous_truth={amb},theirs={round((tp + uth) / N, 4)},strict={round(st / N, 4)}")
        print(line, flush=True); out.write(line + "\n")
