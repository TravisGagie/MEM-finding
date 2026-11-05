#!/usr/bin/env python3
"""Reads Cliffy's Figure 2 (classification time by biome and region) from the vector graphics of its PDF.

    pdftocairo -svg -f PAGE -l PAGE PAPER.pdf page.svg
    cliffy_fig2_values.py page.svg BOTTOM_TICK_Y TOP_TICK_Y

Pages: 12 of the bioRxiv preprint (ticks 425.3 and 323.6) and 9 of the Genome Research paper (ticks 195.02 and 89.74).
The y axis has 9 ticks, 500 s to 4500 s; the preprint's labels are text, and the journal's are outlines, so its ticks are
taken to be the same.  Series come from the figure's own legend: dashed, DNA minimizers; dotted, minimizers; dot-dash, no
digestion; red, Cliffy; teal and solid, Kraken 2.  exp3.r's linetype values don't match these.
Panels, left to right: aquatic, human gut, soil.  Regions, left to right: V1-V2, V3-V4, V4, V4-V5.
"""
import re, sys
s = open(sys.argv[1]).read(); y0 = float(sys.argv[2]); y1 = float(sys.argv[3])
val = lambda y: 500 + (y0 - y) * 4000 / (y0 - y1)
L = []
for p in re.findall(r'<path[^>]*?/>', s):
    if 'fill="none"' not in p or " C " in p: continue
    d = re.search(r' d="([^"]*)"', p).group(1); pts = [(float(x), float(y)) for x, y in re.findall(r'[ML] ([-\d.]+) ([-\d.]+)', d)]
    if len(pts) != 4 or len({q[1] for q in pts}) < 2 and False: continue
    t = re.search(r'transform="matrix\(([^)]*)\)"', p); a, b, c, e, tx, ty = [float(x) for x in t.group(1).split(",")] if t else (1, 0, 0, 1, 0, 0)
    P = [(a * x + c * y + tx, b * x + e * y + ty) for x, y in pts]
    st = re.search(r'stroke="([^"]+)"', p); da = re.search(r'stroke-dasharray="([^"]+)"', p)
    if not st: continue
    col = "Cliffy" if st.group(1).startswith("rgb(97.25") else "Kraken2" if st.group(1).startswith("rgb(0%, 74.9") else None
    if not col or y1 - 5 > min(q[1] for q in P) or max(q[1] for q in P) > y0 + 12: continue
    if abs(P[1][0] - P[0][0] - (P[3][0] - P[2][0])) > .5 or P[1][0] - P[0][0] < 10: continue
    n = len(da.group(1).split()) if da else 0
    pat = {0: "-", 4: "no digestion"}.get(n)
    if n == 2:
        u = [float(x) for x in da.group(1).split()]; pat = "DNA minimizers" if u[0] == u[1] else "minimizers"
    L.append((round(P[0][0]), col, pat, [val(q[1]) for q in P]))
xs = sorted({x[0] for x in L}); B = ["aquatic", "human gut", "soil"]
print("lines", len(L), "panels at x", xs)
assert len(L) == 12 and len(xs) == 3
print("biome      method   digestion        V1-V2   V3-V4      V4   V4-V5")
for x in L:
    print("%-10s %-8s %-14s " % (B[xs.index(x[0])], x[1], x[2]) + " ".join("%7.0f" % v for v in x[3]))
