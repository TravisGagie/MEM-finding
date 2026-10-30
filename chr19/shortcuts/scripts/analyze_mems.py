#!/usr/bin/env python3
"""Merge and summarize the MEM-finding CSVs of run_mem_experiment.sh.

    analyze_mems.py OUTDIR INDEX REVERSED_INDEX

Writes OUTDIR/mem_results.csv (one row per pattern, configuration,
pattern set, algorithm and minimum length) and OUTDIR/mem_summary.txt,
with the size of both indexes for each configuration and, for each
pattern set, algorithm, minimum length and pattern length, the speedup
(total plain time / total time with shortcuts), the time per pattern,
the MEMs per pattern and the fraction of backward steps saved.
"""
import csv
import glob
import os
import re
import sys
from collections import defaultdict


def size_mib(index, config):
    total = os.path.getsize(os.path.join(index, "rlbwt.bin"))
    total += os.path.getsize(os.path.join(index, config + ".shortcuts"))
    return total / 2**20


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    out, index, reversed_index = sys.argv[1:]
    rows = []
    for path in sorted(glob.glob(os.path.join(out, "mems", "*.csv"))):
        m = re.fullmatch(r"(ell\d+_theta\d+)_([\d.]+)_(\w+)\.csv", os.path.basename(path))
        if not m:
            continue
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                rows.append({"config": m[1], "percent": m[2], **row})
    if not rows:
        sys.exit("no results")
    with open(os.path.join(out, "mem_results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    sums = defaultdict(lambda: defaultdict(float))
    for r in rows:
        k = (r["config"], r["percent"], r["algorithm"], int(r["min_length"]), int(r["length"]))
        s = sums[k]
        s["n"] += 1
        for col in ("plain_ns", "shortcut_ns", "plain_steps", "shortcut_steps", "mems", "shortcuts"):
            s[col] += float(r[col])

    plain_mib = os.path.getsize(os.path.join(index, "rlbwt.bin")) / 2**20
    plain_mib += os.path.getsize(os.path.join(reversed_index, "rlbwt.bin")) / 2**20
    lines = [f"both RLBWTs: {plain_mib:.1f} MiB"]
    configs = sorted({k[0] for k in sums}, key=lambda c: int(re.match(r"ell(\d+)", c)[1]))
    for config in configs:
        total = size_mib(index, config) + size_mib(reversed_index, config)
        lines.append(f"\n{config}: both indexes with shortcuts {total:.1f} MiB (+{100 * (total / plain_mib - 1):.0f}%)")
        order = lambda k: (float(k[1]), k[2] != "fb", k[3], k[4])  # forward-backward first
        for key in sorted((k for k in sums if k[0] == config), key=order):
            _, p, alg, L, length = key
            s = sums[key]
            lines.append(
                f"  {p}% substituted, {alg}, L = {L}, length {length}: speedup {s['plain_ns'] / s['shortcut_ns']:.2f}, "
                f"{s['plain_ns'] / s['n'] / 1e3:.1f} -> {s['shortcut_ns'] / s['n'] / 1e3:.1f} us per pattern, "
                f"{s['mems'] / s['n']:.1f} MEMs, {s['shortcuts'] / s['n']:.1f} shortcuts, "
                f"{100 * (1 - s['shortcut_steps'] / s['plain_steps']):.0f}% of backward steps saved")
    with open(os.path.join(out, "mem_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
