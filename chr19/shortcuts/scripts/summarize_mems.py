#!/usr/bin/env python3
"""Summarize query_mems CSVs, one per configuration.

    summarize_mems.py CSV...

For each configuration (the file name, and the config column if any),
algorithm and minimum length, prints the speedup (total
plain time / total time with shortcuts), the backward steps per pattern
character with and without shortcuts, the shortcuts per 1000 characters,
and the number and mean length of the MEMs found.
"""
import csv
import os
import sys
from collections import defaultdict


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    print("config\talgorithm\tL\tpatterns\tplain_s\tshortcut_s\tspeedup\tplain_steps_per_char\t"
          "shortcut_steps_per_char\tshortcuts_per_kchar\tMEMs\tmean_MEM_length")
    for path in sys.argv[1:]:
        sums = defaultdict(lambda: defaultdict(float))
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                stem = os.path.basename(path).rsplit(".", 1)[0]
                config = stem + ("/" + r["config"] if r.get("config") else "")
                s = sums[(config, r["algorithm"], int(r["min_length"]))]
                s["n"] += 1
                for col in ("length", "mems", "mem_chars", "plain_ns", "shortcut_ns", "plain_steps",
                            "shortcut_steps", "shortcuts"):
                    s[col] += float(r[col])
        for (config, alg, L), s in sorted(sums.items()):
            print(f"{config}\t{alg}\t{L}\t{int(s['n'])}\t{s['plain_ns'] / 1e9:.1f}\t{s['shortcut_ns'] / 1e9:.2f}\t"
                  f"{s['plain_ns'] / s['shortcut_ns']:.2f}\t{s['plain_steps'] / s['length']:.3f}\t"
                  f"{s['shortcut_steps'] / s['length']:.3f}\t{1000 * s['shortcuts'] / s['length']:.2f}\t"
                  f"{int(s['mems'])}\t{s['mem_chars'] / max(1, s['mems']):.0f}")


if __name__ == "__main__":
    main()
