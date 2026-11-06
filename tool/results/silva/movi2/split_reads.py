#!/usr/bin/env python3
"""Reads to Movi 2: FASTA, split at non-ACGT characters (as katka find does), keeping segments of at least L, named name:offset."""
import re, sys
L = int(sys.argv[3]); o = open(sys.argv[2], "w")
for i, l in enumerate(open(sys.argv[1])):
    if i % 4 == 0: name = l[1:].split()[0]
    elif i % 4 == 1:
        for m in re.finditer(r"[ACGT]+", l.strip().upper()):
            if m.end() - m.start() >= L: o.write(">%s:%d\n%s\n" % (name, m.start(), m.group()))
