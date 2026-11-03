#!/usr/bin/env python3
# brute_force.py DUMP J D M ...: recomputes rlz_parse's counts naively from its phrase dump (small data only):
# H0, H1, greedy RLZ, GDC-style RLZ and its relative-pointer runs, and RLZAP-style parses, with the same rules
import sys, math, collections
L = open(sys.argv[1]).read().split('\n'); n, R, sigma = map(int, L[0].split())
ids = [int(l.split()[1]) for l in L[1:] if l]; m = len(ids); T = ids[:R]
c = collections.Counter(ids); h0 = -sum(v * math.log2(v / m) for v in c.values()) / m
pc = collections.Counter(zip(ids, ids[1:])); cc = collections.Counter(ids[:-1])
h1 = -sum(v * math.log2(v / cc[a]) for (a, b), v in pc.items()) / (m - 1)
print('H0 = %.4f, H1 = %.4f' % (h0, h1))
def lce(p, q):
    l = 0
    while p + l < m and q + l < R and ids[p + l] == T[q + l]: l += 1
    return l
first = collections.defaultdict(list)
for q in range(R): first[T[q]].append(q)
def longest(p):        # length and the sorted set of sources achieving it (every q with T[q] = ids[p])
    best, srcs = 0, []
    for q in first.get(ids[p], []):
        l = lce(p, q)
        if l > best: best, srcs = l, [q]
        elif l == best and l > 0: srcs.append(q)
    return best, srcs
z = lit = 0; p = R
while p < m:
    d, _ = longest(p); z += 1
    if d == 0: lit += 1; p += 1
    else: p += d
print('RLZ: %d factors, %d literals' % (z, lit))
zg = runs = 0; prev = None; p = R
while p < m:
    d, srcs = longest(p)
    if d > 0:
        q = min(srcs, key=lambda x: (x != (p + prev if prev is not None else -1), x))     # continue prev if possible
        # rlz_parse takes SA[lo] (the lexicographically smallest suffix) when it cannot continue; compare only runs/factors
        rel = q - p
        if rel != prev: runs += 1; prev = rel
    zg += 1; p += d + 1
print('GDC: %d factors, %d runs' % (zg, runs))
for a in range(2, len(sys.argv), 3):
    J, D, M = int(sys.argv[a]), int(sys.argv[a + 1]), int(sys.argv[a + 2])
    zf = lits = changes = expl = 0; have = False; rel = 0; p = R
    while p < m:
        done = False
        if have:
            for j in range(J + 1):
                if done or p + j >= m: break
                for dl in [0] + [x for k in range(1, D + 1) for x in (k, -k)]:
                    src = p + j + rel + dl
                    if 0 <= src < R and lce(p + j, src) >= M:
                        zf += 1; lits += j
                        if dl: changes += 1; rel += dl
                        p += j + lce(p + j, src); done = True; break
        if done: continue
        d, srcs = longest(p); zf += 1
        if d == 0: lits += 1; p += 1; continue
        nr = (p + rel) - p if have and (p + rel) in srcs else None
        if nr is None: nr = min(srcs) - p      # stand-in for SA[lo]; affects only which pointer, so counts may differ slightly
        if not have or nr != rel: expl += 1; rel = nr; have = True
        p += d
    print('RLZAP J=%d D=%d M=%d: %d factors, %d literals, %d changes, %d explicit' % (J, D, M, zf, lits, changes, expl))
