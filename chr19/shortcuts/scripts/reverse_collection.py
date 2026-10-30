#!/usr/bin/env python3
"""Write the reversed collection: the lines of a collection file in reverse
order, each reversed, as `tac IN | rev` would (text.hpp).

    reverse_collection.py IN OUT

IN must have one sequence per line, with no blank lines and a final
newline, as prepare_collection.py writes. Its BWT, from Big-BWT, can be
given to build_index --reversed --bwt.
"""
import mmap
import os
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    with open(sys.argv[1], "rb") as f, open(sys.argv[2], "wb") as out:
        size = os.fstat(f.fileno()).st_size
        if size == 0:
            return
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        if mm[size - 1:size] != b"\n":
            sys.exit("the input must end with a newline")
        end = size - 1  # the newline ending the current line
        while end >= 0:
            start = mm.rfind(b"\n", 0, end) + 1
            if start == end:
                sys.exit("the input must not have blank lines")
            out.write(mm[start:end][::-1])
            out.write(b"\n")
            end = start - 1


if __name__ == "__main__":
    main()
