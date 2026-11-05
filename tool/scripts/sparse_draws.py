#!/usr/bin/env python3
"""The sparse dial's abundance tables: the aquatic sample's 100 abundances on 100 random genera (tool/results/silva/README.md,
"The sparse dial").

    sparse_draws.py GENERA.tsv AQUATIC.tsv OUT [--draws 3] [--seed 1]

GENERA.tsv is scripts/sparse_genera.py's.  A genus is usable in a region if its name is unique and not "uncultured", it
has the region, and over 55% of its regions are ones ART can read (scripts/sparse_genera.py's nlong; MicrobeMixer's ART
coverage is twice what it needs).  Each region gets its own draw: only about a third of SILVA's genera have V1-V2, so one draw
across the regions would leave most genera empty in some region, or, restricted to the genera usable in all four, would
favour large genera (median N_g 21, against 4).  The k-th largest abundance goes to the k-th genus drawn.

Writes OUT/uniform.<k>/<region>.tsv (all usable genera) and OUT/small.<k>/<region>.tsv (the usable genera with N_g <= 2,
the smallest third of SILVA's genera, ties included: 3,974 genera against 3,039), and OUT/draws.txt.
"""
import argparse, collections, os, random, statistics
ap = argparse.ArgumentParser()
ap.add_argument("genera"); ap.add_argument("aquatic"); ap.add_argument("out")
ap.add_argument("--draws", type=int, default=3); ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
REGIONS = "V1_V2 V3_V4 V4_V4 V4_V5".split()
rows = [l.rstrip("\n").split("\t") for l in open(a.genera)][1:]
names = [r[2].rstrip(";").split(";")[-1] for r in rows]; cnt = collections.Counter(names)
ab = list(open(a.aquatic)); head, ab = ab[0], [l.rstrip("\n").split("\t") for l in ab[1:]]
assert len(ab) == 100
rnd = random.Random(a.seed)
log = open(os.path.join(a.out, "draws.txt"), "w") if os.path.isdir(a.out) else None
if log is None:
    os.makedirs(a.out); log = open(os.path.join(a.out, "draws.txt"), "w")
for name, cap in (("uniform", None), ("small", 2)):
    for j, r in enumerate(REGIONS):
        use = [x for x, n in zip(rows, names) if cnt[n] == 1 and n != "uncultured" and int(x[3 + 2 * j]) > 0
               and int(x[4 + 2 * j]) / int(x[3 + 2 * j]) > 0.55 and (cap is None or int(x[1]) <= cap)]
        print("%s %s: %d usable genera, median N_g %g" % (name, r, len(use), statistics.median(int(x[1]) for x in use)), file=log)
        assert len(use) >= 100, (name, r, len(use))
    for k in range(1, a.draws + 1):
        d = os.path.join(a.out, "%s.%d" % (name, k)); os.makedirs(d, exist_ok=True)
        for j, r in enumerate(REGIONS):
            use = [x for x, n in zip(rows, names) if cnt[n] == 1 and n != "uncultured" and int(x[3 + 2 * j]) > 0
                   and int(x[4 + 2 * j]) / int(x[3 + 2 * j]) > 0.55 and (cap is None or int(x[1]) <= cap)]
            pick = rnd.sample(use, 100)
            print("%s.%d %s: N_g median %g, range %d to %d" % (name, k, r, statistics.median(int(x[1]) for x in pick),
                  min(int(x[1]) for x in pick), max(int(x[1]) for x in pick)), file=log)
            with open(os.path.join(d, r + ".tsv"), "w") as o:
                o.write(head)
                for g, l in zip(pick, ab): o.write("%s\t%s\t%s\n" % (g[2], l[1], l[2]))
log.close()
