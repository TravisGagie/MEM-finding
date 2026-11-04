#!/usr/bin/env python3
"""The tags of the tag runs, as 32-bit integers, for bigrepair -i (tool/NOTES.md, "Grammar tags").

    tagruns_to_int32.py PREFIX.tagruns OUT.int32

PREFIX.tagruns (or PREFIX.sp.tagruns with --by-tag) is the cache `katka index --tags` writes: "KTRUNS01", n
(u64), the run starts (u64 count, then u64 values) and the runs' tags (u64 count, then u32 values).  OUT gets
the tags in run order, one little-endian 32-bit integer per run: the sequence the tool's built-in RePair
compresses, so `katka index --tags --bigrepair OUT` then reads bigrepair's OUT.R and OUT.C.  Streams the file
(16 MB at a time), so SILVA's 359 million runs need no memory to speak of.
"""
import struct, sys

if len(sys.argv) != 3:
    sys.exit(__doc__)
with open(sys.argv[1], "rb") as f, open(sys.argv[2], "wb") as o:
    if f.read(8) != b"KTRUNS01":
        sys.exit("%s is not a katka tag-run cache (KTRUNS01)" % sys.argv[1])
    (n,) = struct.unpack("<Q", f.read(8))
    (k,) = struct.unpack("<Q", f.read(8))
    f.seek(8 * k, 1)                                   # the run starts
    (k2,) = struct.unpack("<Q", f.read(8))
    if k2 != k:
        sys.exit("bad cache: %d starts but %d tags" % (k, k2))
    left, mx = k, 0
    while left:
        b = f.read(4 * min(left, 1 << 22))
        if not b:
            sys.exit("short read")
        vals = struct.unpack("<%dI" % (len(b) // 4), b)
        mx = max(mx, max(vals))
        if mx >= 1 << 30:
            sys.exit("a tag is >= 2^30, which bigrepair -i does not take")
        o.write(b)
        left -= len(b) // 4
print("%s: %d runs over n = %d rows, tags 0..%d" % (sys.argv[2], k, n, mx), file=sys.stderr)
