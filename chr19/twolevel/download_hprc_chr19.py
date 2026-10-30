#!/usr/bin/env python3
"""Download the chromosome 19 dataset of "Reportable Shortcuts with RLBWTs".

    python3 download_hprc_chr19.py OUTDIR [--lines] [--raw]

The dataset comes from release 2 of the Human Pangenome Reference
Consortium (HPRC): chromosome 19 of CHM13 (the reference), of both
haplotypes of 16 individuals (the collection), and of haplotype 1 of two
further individuals, HG00423 and HG00438 (the queries). The 16 individuals
are the first by sample ID whose two haplotypes each have chromosome 19 as
a single sequence of at least 90% of the length of CHM13's, leaving out
the two query individuals. All sequences are in the orientation of CHM13.

Only the bases of chromosome 19 are downloaded (about 1 GB in all), block
by block from the bgzipped assemblies on S3, using their .fai and .gzi
indexes; the locations of the assemblies come from the HPRC release 2
index on GitHub. Needs Python 3.8 or later and nothing else.

Writes OUTDIR/collection.fa (CHM13 first, then the 32 haplotypes) and
OUTDIR/queries.fa. As in the paper, sequences are upper-cased and every
character other than A, C, G and T becomes N; with --raw, they are left
as in the assemblies. With --lines, also writes OUTDIR/collection.txt and
OUTDIR/queries.txt, one sequence per line, which is the input format of
the paper's code. Each sequence is checked against its length and number
of Ns in the paper's dataset. Rerunning skips sequences already written.
"""
import argparse
import csv
import io
import os
import struct
import sys
import time
import urllib.request
import zlib

INDEX_URL = ("https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/"
             "main/data_tables/assemblies_release2_v1.0.index.csv")
S3_HTTPS = "https://s3-us-west-2.amazonaws.com/"

# (role, sample, haplotype, sequence name in the assembly, length, number of N
# after converting characters other than A/C/G/T)
DATASET = [
    ("reference", "CHM13", "0", "CHM13#0#chr19", 61707364, 0),
    ("collection", "HG00097", "1", "HG00097#1#CM094070.1", 61077493, 0),
    ("collection", "HG00097", "2", "HG00097#2#CM094085.1", 61125280, 142121),
    ("collection", "HG00099", "1", "HG00099#1#CM087326.1", 62336960, 0),
    ("collection", "HG00099", "2", "HG00099#2#CM087370.1", 62320705, 0),
    ("collection", "HG00126", "1", "HG00126#1#CM090118.1", 61924483, 20782),
    ("collection", "HG00126", "2", "HG00126#2#CM090133.1", 61592150, 0),
    ("collection", "HG00128", "1", "HG00128#1#CM090090.1", 60445756, 18644),
    ("collection", "HG00128", "2", "HG00128#2#CM090105.1", 61467610, 0),
    ("collection", "HG00140", "1", "HG00140#1#CM087120.1", 61951028, 0),
    ("collection", "HG00140", "2", "HG00140#2#CM087145.1", 59211742, 20814),
    ("collection", "HG00232", "1", "HG00232#1#CM090000.1", 60318481, 0),
    ("collection", "HG00232", "2", "HG00232#2#CM090032.1", 61648317, 0),
    ("collection", "HG00235", "1", "HG00235#1#CM094378.1", 61759272, 102738),
    ("collection", "HG00235", "2", "HG00235#2#CM094405.1", 61339476, 0),
    ("collection", "HG00272", "1", "HG00272#1#CM094199.1", 61221721, 0),
    ("collection", "HG00272", "2", "HG00272#2#CM094225.1", 63486065, 100),
    ("collection", "HG00280", "1", "HG00280#1#CM087157.1", 62168755, 0),
    ("collection", "HG00280", "2", "HG00280#2#CM087169.1", 59604861, 0),
    ("collection", "HG00290", "1", "HG00290#1#CM089990.1", 61006069, 100),
    ("collection", "HG00290", "2", "HG00290#2#CM090011.1", 61094987, 0),
    ("collection", "HG00320", "1", "HG00320#1#CM089959.1", 61987789, 100),
    ("collection", "HG00320", "2", "HG00320#2#CM089975.1", 61506238, 260601),
    ("collection", "HG00321", "1", "HG00321#1#CM089551.1", 60983011, 0),
    ("collection", "HG00321", "2", "HG00321#2#CM089581.1", 63745078, 0),
    ("collection", "HG00329", "1", "HG00329#1#CM094236.1", 60473049, 0),
    ("collection", "HG00329", "2", "HG00329#2#CM094247.1", 62456877, 0),
    ("collection", "HG00344", "1", "HG00344#1#CM094509.1", 62001390, 419721),
    ("collection", "HG00344", "2", "HG00344#2#CM094520.1", 60540825, 0),
    ("collection", "HG00350", "1", "HG00350#1#CM094479.1", 61981818, 0),
    ("collection", "HG00350", "2", "HG00350#2#CM094499.1", 64240453, 34866),
    ("collection", "HG00408", "1", "HG00408#1#CM085963.1", 61061584, 0),
    ("collection", "HG00408", "2", "HG00408#2#CM085976.1", 62956045, 0),
    ("query", "HG00423", "1", "HG00423#1#CM089409.1", 61651182, 163509),
    ("query", "HG00438", "1", "HG00438#1#CM089182.1", 60966508, 19557),
]


def _acgt_table():
    table = bytearray(b"N" * 256)
    for c in b"ACGT":
        table[c] = c
        table[c + 32] = c
    return bytes(table)


ACGT = _acgt_table()


def https(url):
    return S3_HTTPS + url[len("s3://"):] if url.startswith("s3://") else url


def fetch(url, start=None, end=None, tries=5):
    """The bytes of url (or a local path), or of its byte range [start, end)."""
    if "://" not in url:
        with open(url, "rb") as f:
            f.seek(start or 0)
            return f.read() if end is None else f.read(end - (start or 0))
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url)
            if start is not None:
                req.add_header("Range", f"bytes={start}-" + ("" if end is None else str(end - 1)))
            with urllib.request.urlopen(req, timeout=120) as resp:
                data = resp.read()
                if start is not None and resp.status == 200:  # the server ignored the range
                    data = data[start:end]
                return data
        except Exception as e:
            if attempt + 1 == tries:
                raise
            print(f"  retrying {url}: {e}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))


def extract(row, name):
    """Sequence `name` of the bgzipped assembly in `row`, downloading only its blocks."""
    fai = {}
    for line in fetch(https(row["assembly_fai"])).decode().splitlines():
        f = line.split("\t")
        fai[f[0]] = tuple(int(x) for x in f[1:5])
    if name not in fai:
        raise RuntimeError(f"{name} is not in {row['assembly']}")
    length, offset, linebases, linewidth = fai[name]
    data = fetch(https(row["assembly_gzi"]))
    (k,) = struct.unpack_from("<Q", data, 0)
    gzi = [(0, 0)] + [struct.unpack_from("<QQ", data, 8 + 16 * i) for i in range(k)]
    u0 = offset
    u1 = offset + (length // linebases) * linewidth + length % linebases
    first = max(i for i in range(len(gzi)) if gzi[i][1] <= u0)
    later = [i for i in range(len(gzi)) if gzi[i][1] >= u1]
    c0, base = gzi[first]
    c1 = gzi[later[0]][0] if later else None
    comp = fetch(https(row["assembly"]), c0, c1)
    out = io.BytesIO()
    while comp:
        d = zlib.decompressobj(31)
        out.write(d.decompress(comp))
        comp = d.unused_data
    seq = out.getvalue()[u0 - base:u1 - base].replace(b"\n", b"").replace(b"\r", b"")
    if len(seq) != length:
        raise RuntimeError(f"{name}: extracted {len(seq)} bases, expected {length}")
    return seq


def write_fasta(f, name, seq, width=80):
    f.write(b">" + name.encode() + b"\n")
    for i in range(0, len(seq), width):
        f.write(seq[i:i + width] + b"\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--lines", action="store_true", help="also write one sequence per line")
    ap.add_argument("--raw", action="store_true", help="keep the sequences as in the assemblies")
    ap.add_argument("--index", default=INDEX_URL, help=argparse.SUPPRESS)  # for testing
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = {}
    for r in csv.DictReader(io.StringIO(fetch(args.index).decode())):
        rows[(r["sample_id"], r["haplotype"])] = r

    # Each sequence goes to its own file first, so that a rerun can skip it.
    parts = os.path.join(args.outdir, "parts")
    os.makedirs(parts, exist_ok=True)
    for k, (role, sample, hap, name, length, ns) in enumerate(DATASET):
        path = os.path.join(parts, name.replace("#", "_") + ".seq")
        if os.path.exists(path) and os.path.getsize(path) == length:
            continue
        print(f"[{k + 1}/{len(DATASET)}] {name} ({role})", file=sys.stderr)
        seq = extract(rows[(sample, hap)], name)
        if not args.raw:
            seq = seq.translate(ACGT)
            if seq.count(b"N") != ns:
                raise RuntimeError(f"{name}: {seq.count(b'N')} Ns, expected {ns}")
        with open(path + ".tmp", "wb") as f:
            f.write(seq)
        os.replace(path + ".tmp", path)

    outputs = [("collection", ("reference", "collection")), ("queries", ("query",))]
    for stem, roles in outputs:
        entries = [d for d in DATASET if d[0] in roles]
        with open(os.path.join(args.outdir, stem + ".fa"), "wb") as fa:
            for role, sample, hap, name, length, ns in entries:
                with open(os.path.join(parts, name.replace("#", "_") + ".seq"), "rb") as f:
                    write_fasta(fa, name, f.read())
        if args.lines:
            with open(os.path.join(args.outdir, stem + ".txt"), "wb") as txt:
                for role, sample, hap, name, length, ns in entries:
                    with open(os.path.join(parts, name.replace("#", "_") + ".seq"), "rb") as f:
                        txt.write(f.read() + b"\n")
    for role, sample, hap, name, length, ns in DATASET:
        os.remove(os.path.join(parts, name.replace("#", "_") + ".seq"))
    os.rmdir(parts)
    n_coll = sum(1 for d in DATASET if d[0] != "query")
    print(f"wrote {args.outdir}/collection.fa ({n_coll} sequences) and {args.outdir}/queries.fa "
          f"({len(DATASET) - n_coll} sequences)", file=sys.stderr)


if __name__ == "__main__":
    main()
