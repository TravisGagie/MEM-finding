#!/usr/bin/env python3
"""Reads Cliffy's Figure 3 (aquatic, accuracy by clade) from the vector graphics of its PDF.

    pdftocairo -svg -f PAGE -l PAGE PAPER.pdf page.svg
    cliffy_fig3_values.py page.svg

Pages: 9 of the Genome Research paper, 13 of the bioRxiv preprint.  Each panel's y axis is calibrated from its 6 tick marks, 1.00 at the top to 0.75 (the labels of the preprint;
in the journal's PDF they are outlines, so its ticks are assumed the same, and the two versions' values should agree).  Cliffy's series (the red lines of exp3.r's plot) are told apart by dash (solid: LCA;
dashed: approximate document listing) and by the marker at each vertex (triangle: no digestion; square: DNA
minimizers; diamond: minimizers).  Prints the value at each clade.
"""
import re, sys
svg = open(sys.argv[1]).read()
def parse_d(d):
    return [(float(x), float(y)) for x, y in re.findall(r'[ML] ([-\d.]+) ([-\d.]+)', d)]
paths = []
for p in re.findall(r'<path[^>]*?/>', svg):
    d = re.search(r' d="([^"]*)"', p).group(1)
    t = re.search(r'transform="matrix\(([^)]*)\)"', p)
    a, b, c, e, tx, ty = [float(x) for x in t.group(1).split(",")] if t else (1, 0, 0, 1, 0, 0)
    pts = [(a * x + c * y + tx, b * x + e * y + ty) for x, y in parse_d(d)]
    st = re.search(r'stroke="([^"]+)"', p); fl = re.search(r'fill="([^"]+)"', p); da = re.search(r'stroke-dasharray="([^"]+)"', p)
    paths.append(dict(pts=pts, stroke=st.group(1) if st else None, fill=fl.group(1) if fl else None, dash=da.group(1) if da else None, curves=" C " in d))
RED = lambda s: s and s.startswith("rgb(97.25")
lines = [p for p in paths if RED(p["stroke"]) and p["fill"] == "none" and len(p["pts"]) == 6]
markers = [p for p in paths if RED(p["stroke"]) and len(p["pts"]) in (3, 4, 5)]
tm = [p for p in paths if p["stroke"] and p["stroke"].startswith("rgb(19.99") and len(p["pts"]) == 2]
htm = [p for p in tm if abs(p["pts"][0][1] - p["pts"][1][1]) < 1e-6 and 1 < abs(p["pts"][0][0] - p["pts"][1][0]) < 5]
# the y axes: each is 6 ticks, 1.00 (top) to 0.75, as the preprint's labels show; the journal's labels are outlines
cols = {}
for p in htm: cols.setdefault(round(min(q[0] for q in p["pts"]), 0), []).append(round(p["pts"][0][1], 2))
axes = []
for ys in cols.values():
    ys.sort()
    axes += [ys[i:i + 6] for i in range(0, len(ys) - 5, 6) if abs((ys[i + 5] - ys[i]) - 5 * (ys[i + 1] - ys[i])) < .05]
axes = [v for v in axes if any(min(q[1] for q in l["pts"]) >= v[0] - 30 and max(q[1] for q in l["pts"]) <= v[5] + 30 for l in lines)]
xs = sorted({round(l["pts"][0][0], 0) for l in lines})
print("first-vertex x", xs, "; lines", len(lines), "; markers", len(markers), "; axes", [(v[0], v[5]) for v in axes])
assert len(lines) == 24
NAMES = ["V1_V2", "V3_V4", "V4_V4", "V4_V5"]
def panel(l):
    x = xs.index(round(l["pts"][0][0], 0))
    if len(xs) == 4: return NAMES[x], axes[0]
    y = 0 if min(q[1] for q in l["pts"]) < axes[1][0] - 15 else 1
    return NAMES[2 * y + x], axes[y]
def value(y, ax):
    return 1.00 - 0.25 * (y - ax[0]) / (ax[5] - ax[0])
out = {}
for l in lines:
    name, ax = panel(l); v0 = l["pts"][1]              # the marker is matched at Family, away from other series' overlaps at Genus
    best = None
    for m in markers:
        u = sorted(set((round(q[0], 2), round(q[1], 2)) for q in m["pts"]))
        cx = sum(q[0] for q in u) / len(u); cy = sum(q[1] for q in u) / len(u)
        if abs(cx - v0[0]) < 2.5 and abs(cy - v0[1]) < 2.5 and len(u) in (3, 4):
            d = (cx - v0[0]) ** 2 + (cy - v0[1]) ** 2
            if best is None or d < best[0]: best = (d, u)
    shape = "?"
    if best:
        u = best[1]
        shape = "none" if len(u) == 3 else ("dna_minimizers" if len({q[0] for q in u}) == 2 else "minimizers")
    k = (name, "listing" if l["dash"] else "LCA", shape)
    assert k not in out, k
    out[k] = [value(q[1], ax) for q in l["pts"]]
print("\nregion  query    digestion         Genus  Family   Order   Class  Phylum  Domain")
for k in sorted(out):
    print("%-7s %-8s %-15s " % k + " ".join("%7.2f" % (100 * v) for v in out[k]))
print("\nGenus, by query and digestion: V1_V2 V3_V4 V4_V4 V4_V5")
for q in ("LCA", "listing"):
    for sh in ("none", "minimizers", "dna_minimizers"):
        v = [out[(r, q, sh)][0] * 100 for r in NAMES]
        print("  %-8s %-15s %s   range %.1f to %.1f" % (q, sh, " ".join("%.2f" % x for x in v), min(v), max(v)))
