#!/usr/bin/env python3
"""Summarize OUTDIR/results.csv and OUTDIR/index/build.csv.

    analyze.py OUTDIR

Writes OUTDIR/summary.txt with, for each configuration (ell, theta), the
index size and, for each pattern set and length, the speedup (total
plain time / total time with shortcuts, both from the same run), the
shortcuts per pattern and the fraction of matched characters they
cross. If matplotlib is available, also writes OUTDIR/pareto.pdf: index
size against total search time for exact patterns, one panel per
pattern length.
"""
import csv
import os
import sys
from collections import defaultdict


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out = sys.argv[1]
    with open(os.path.join(out, "index", "build.csv"), newline="") as f:
        builds = {(int(r["ell"]), int(r["theta"])): r for r in csv.DictReader(f)}

    # (ell, theta, percent, length) -> sums over the patterns
    sums = defaultdict(lambda: defaultdict(float))
    with open(os.path.join(out, "results.csv"), newline="") as f:
        for r in csv.DictReader(f):
            s = sums[(int(r["ell"]), int(r["theta"]), r["percent"], int(r["length"]))]
            s["patterns"] += 1
            for col in ("plain_ns", "shortcut_ns", "shortcuts", "chars_crossed", "matched", "found"):
                s[col] += float(r[col])

    configs = sorted({k[:2] for k in sums})
    percents = sorted({k[2] for k in sums}, key=float)
    lengths = sorted({k[3] for k in sums})
    any_build = next(iter(builds.values()))
    rlbwt_mib = int(any_build["rlbwt_bytes"]) / 2**20

    lines = [f"RLBWT: n = {any_build['n']}, r = {any_build['runs']}, {rlbwt_mib:.1f} MiB"]
    for ell, theta in configs:
        b = builds[(ell, theta)]
        total = int(b["total_bytes"]) / 2**20
        lines.append(f"\nell = {ell}, theta = {theta}: {total:.1f} MiB (+{100 * (total / rlbwt_mib - 1):.0f}%), "
                     f"{b['entrances']} entrances")
        for p in percents:
            for L in lengths:
                s = sums[(ell, theta, p, L)]
                lines.append(f"  {p}% substituted, length {L}: speedup {s['plain_ns'] / s['shortcut_ns']:.2f}, "
                             f"{s['shortcuts'] / s['patterns']:.1f} shortcuts per pattern, "
                             f"{100 * s['chars_crossed'] / max(1, s['matched']):.0f}% of matched characters crossed, "
                             f"{int(s['found'])} of {int(s['patterns'])} found")
    with open(os.path.join(out, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    exact = percents[0]
    fig, axes = plt.subplots(1, len(lengths), figsize=(2.4 * len(lengths), 2.6))
    for ax, L in zip(axes, lengths):
        for ell in sorted({c[0] for c in configs}):
            pts = [(int(builds[(e, t)]["total_bytes"]) / 2**20, sums[(e, t, exact, L)]["shortcut_ns"] / 1e6)
                   for e, t in configs if e == ell]
            ax.scatter(*zip(*pts), s=10, label=f"$\\ell={ell}$")
        plain_ms = sorted(sums[(e, t, exact, L)]["plain_ns"] / 1e6 for e, t in configs)[len(configs) // 2]
        ax.scatter([rlbwt_mib], [plain_ms], marker="*", s=60, color="black", label="RLBWT only")
        ax.set_title(f"length {L}", fontsize=8)
        ax.set_xlabel("index size (MiB)", fontsize=7)
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=6)
    axes[0].set_ylabel("query time (ms)", fontsize=7)
    axes[0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "pareto.pdf"))


if __name__ == "__main__":
    main()
