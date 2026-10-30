#!/usr/bin/env python3
"""Merge the per-configuration query CSVs into OUTDIR/results.csv.

    collect_results.py OUTDIR

Each row of results.csv is one search: the configuration (ell, theta),
the substitution percentage of the pattern set, and query_index's
columns.
"""
import csv
import glob
import os
import re
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out = sys.argv[1]
    writer = None
    with open(os.path.join(out, "results.csv"), "w", newline="") as f:
        for path in sorted(glob.glob(os.path.join(out, "queries", "*.csv"))):
            m = re.fullmatch(r"ell(\d+)_theta(\d+)_([\d.]+)\.csv", os.path.basename(path))
            if not m:
                continue
            with open(path, newline="") as g:
                reader = csv.DictReader(g)
                if writer is None:
                    writer = csv.DictWriter(f, ["ell", "theta", "percent"] + reader.fieldnames)
                    writer.writeheader()
                for row in reader:
                    writer.writerow({"ell": m[1], "theta": m[2], "percent": m[3], **row})


if __name__ == "__main__":
    main()
