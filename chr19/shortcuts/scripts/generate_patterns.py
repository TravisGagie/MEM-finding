#!/usr/bin/env python3
"""Generate the patterns used in the experiments.

    generate_patterns.py COLLECTION OUTDIR [PERCENTS]

Takes 200 substrings of each length 200, 400, ..., 1000 from the
collection file, starting at uniformly random positions and consisting
only of A, C, G and T (so none crosses from one line to the next), and
writes them as OUTDIR/patterns_<p>.txt with each character substituted
with probability p% for each p in PERCENTS (default 0,0.1,1). Line i is
the same substring in every file. The random seed is fixed, so the
patterns depend only on the collection and PERCENTS. The experiments on
the bacterial collection used PERCENTS = 0,0.1,1, and those on cere
PERCENTS = 0,1.
"""
import mmap
import os
import random
import sys

LENGTHS = [200, 400, 600, 800, 1000]
PER_LENGTH = 200
SEED = 20260913


def substitute(s, rate):
    if rate == 0:
        return s.decode("ascii")
    out = bytearray(s)
    for i in range(len(out)):
        if random.random() < rate:
            out[i] = ord(random.choice([c for c in "ACGT" if ord(c) != out[i]]))
    return out.decode("ascii")


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    collection, outdir = sys.argv[1], sys.argv[2]
    percents = (sys.argv[3] if len(sys.argv) == 4 else "0,0.1,1").split(",")
    os.makedirs(outdir, exist_ok=True)
    random.seed(SEED)
    with open(collection, "rb") as f:
        text = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    patterns = {p: [] for p in percents}
    for length in LENGTHS:
        for _ in range(PER_LENGTH):
            while True:
                start = random.randint(0, len(text) - length)
                s = text[start:start + length]
                if not s.translate(None, b"ACGT"):  # only A, C, G, T
                    break
            for p in percents:
                patterns[p].append(substitute(s, float(p) / 100))
    for p, lines in patterns.items():
        with open(os.path.join(outdir, f"patterns_{p}.txt"), "w") as out:
            out.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
