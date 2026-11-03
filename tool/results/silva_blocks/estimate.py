#!/usr/bin/env python3
# Backward steps per read with order-k character-level blocks, against the 486 of katka find now
# (preset silva: L = 30, lookup table t = 10, RLCSA), on the same 8,000 reads. Two estimates:
#  A. from katka-blocks (d), which classifies every step of its own BML run (no lookup table: 677 steps per
#     read): remove the steps inside long crossings, add the block steps replacing them, and start from 486
#     (the table saves the first 10 steps of each search, at its right end, which the model charges as
#     character steps anyway).
#  B. from the bounds of ../silva_profile/step_bounds.py (split at 64, to match katka-blocks' length classes):
#     crossing steps of MEMs <= 64 and >= 65 scaled by katka-blocks' (c) steps per character for those
#     lengths (character-weighted over its classes); failed windows unchanged.
import re, sys
R = 8000; NOW = 486.1; TOTAL_NOW = 3889139
lens = [int(l.split('\t')[2]) - int(l.split('\t')[1]) for l in open(sys.argv[1])]    # profile/mems.out
chars = {'17-32': 0, '33-64': 0, '65-128': 0, '129-256': 0}
for x in lens: chars['17-32' if x <= 32 else '33-64' if x <= 64 else '65-128' if x <= 128 else '129-256'] += x
bounds = {'short': (1699391, 2087976), 'long': (1333243, 2163754)}         # profile/step_bounds64.txt
print('k   A: steps/read (x fewer)   B: steps/read (x fewer)   rho <=64  rho >=65')
for k, f in ((8, sys.argv[2]), (12, sys.argv[2]), (16, sys.argv[3])):
    t = open(f).read(); sec = t[t.index('== k = %d\n' % k):]; sec = sec[:sec.index('total under the model')]
    rho = {c: float(re.search(r'MEM length %s\s+\d+ MEMs, ([0-9.]+) steps' % c, sec).group(1)) for c in chars}
    inside = int(re.search(r'inside long crossings:\s+(\d+)', sec).group(1)); blocks = int(re.search(r'-> (\d+) block steps', sec).group(1))
    A = NOW - inside / R + blocks / R
    rs = (rho['17-32'] * chars['17-32'] + rho['33-64'] * chars['33-64']) / (chars['17-32'] + chars['33-64'])
    rl = (rho['65-128'] * chars['65-128'] + rho['129-256'] * chars['129-256']) / (chars['65-128'] + chars['129-256'])
    # B: failed = total - crossing; with blocks, total' = failed + rs*short + rl*long = total - (1-rs)*short - (1-rl)*long
    lo = (TOTAL_NOW - (1 - rs) * bounds['short'][1] - (1 - rl) * bounds['long'][1]) / R
    hi = (TOTAL_NOW - (1 - rs) * bounds['short'][0] - (1 - rl) * bounds['long'][0]) / R
    print('%2d  %6.1f (%.2fx)            %5.1f-%5.1f (%.2f-%.2fx)       %.3f     %.3f' % (k, A, NOW / A, lo, hi, NOW / hi, NOW / lo, rs, rl))
