#!/usr/bin/env python3
"""Turn FASTA files into a collection file with the reference first.

build_index takes the reference R as a prefix of the collection
(--reference-length N). For a collection with
several species, a single prefix would cover only one species, so this
script writes the representative genomes first and everything else after
them. R is then the representatives, back to back; the record separators
between them are not A/C/G/T, so no window or tier string crosses a
boundary between two representatives.

Output: one record per line (no headers), upper-cased, with every
character other than A, C, G, T replaced by N. build_index reads such a
file as one sequence per line and puts a separator after each.

Usage:
  prepare_collection.py OUT.txt --reps REP1.fa [REP2.fa ...] --rest A.fa [B.fa ...]
  (inputs may be gzipped; a file listed under --reps is not repeated
  from --rest)

  prepare_collection.py OUT.txt --ref-min-len L [--max-bases-per-file B] --rest A.fa [B.fa ...]
  (from EACH --rest file, the first record of length at least L goes into
  the reference -- with one file per species and L = 2000000, that is a
  complete chromosome of each species -- and is not repeated later.
  --max-bases-per-file stops taking records from a file once B bases have
  been taken from it, the reference record included, so that one species
  does not dominate and the collection fits on disk)

Prints the value to pass as --reference-length: the total length of the
representatives' records plus the separators between them.
"""
import gzip
import os
import sys

# Upper-case A/C/G/T; every other byte becomes N. One translate per record.
_tab = bytearray(b"N" * 256)
for _c in b"ACGT":
    _tab[_c] = _c
    _tab[_c + 32] = _c
TABLE = bytes(_tab)


def records(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rb") as f:
        seq = []
        for line in f:
            line = line.rstrip(b"\r\n")
            if line.startswith(b">"):
                if seq:
                    yield b"".join(seq)
                seq = []
            elif line:
                seq.append(line)
        if seq:
            yield b"".join(seq)


def clean(s):
    return s.translate(TABLE)


def main():
    args = sys.argv[1:]
    if not args or ("--reps" not in args and "--ref-min-len" not in args):
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    out = args[0]
    reps, rest, cur = [], [], None
    min_len = None
    max_per_file = None
    it = iter(args[1:])
    for a in it:
        if a == "--ref-min-len":
            min_len = int(next(it))
        elif a == "--max-bases-per-file":
            max_per_file = int(next(it))
        elif a == "--reps":
            cur = reps
        elif a == "--rest":
            cur = rest
        elif cur is None:
            sys.exit(f"unexpected argument {a}")
        else:
            cur.append(a)
    seen = {os.path.realpath(p) for p in reps}
    rest = [p for p in rest if os.path.realpath(p) not in seen]

    ref_len = 0
    n_rep_records = 0
    total = 0
    n_records = 0
    skip = {}  # file index -> index of the record taken from it for the reference
    taken = {}  # file index -> bases taken from it so far
    with open(out, "wb") as w:
        def emit(seq, is_ref):
            nonlocal ref_len, n_rep_records, total, n_records
            w.write(seq + b"\n")
            total += len(seq)
            n_records += 1
            if is_ref:
                ref_len += len(seq)
                n_rep_records += 1

        if min_len is not None:
            for fi, p in enumerate(rest):
                for ri, s in enumerate(records(p)):
                    if len(s) >= min_len:
                        skip[fi] = ri
                        taken[fi] = len(s)
                        emit(clean(s), True)
                        print(f"reference record from {p}: record {ri}, {len(s)} bases", file=sys.stderr)
                        break
                else:
                    sys.exit(f"no record of length at least {min_len} in {p}")
        for p in reps:
            for s in records(p):
                emit(clean(s), True)
        for fi, p in enumerate(rest):
            for ri, s in enumerate(records(p)):
                if skip.get(fi) == ri:
                    continue
                if max_per_file is not None and taken.get(fi, 0) + len(s) > max_per_file:
                    break
                taken[fi] = taken.get(fi, 0) + len(s)
                emit(clean(s), False)
            print(f"{p}: {taken.get(fi, 0)} bases taken", file=sys.stderr)
    ref_prefix = ref_len + max(0, n_rep_records - 1)
    print(f"records: {n_records} ({n_rep_records} in the reference); "
          f"sequence: {total} bases; reference: {ref_len} bases", file=sys.stderr)
    print(ref_prefix)


if __name__ == "__main__":
    main()
