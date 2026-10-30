#!/usr/bin/env python3
"""Extract one chromosome from the HPRC release 2 assemblies.

    prepare_hprc_chr19.py OUTDIR [--chrom chr19] [--individuals 16] [--queries 2]

Reads the release 2 assembly index, takes the chromosome from CHM13 as the
reference, and then goes through the individuals in order of sample ID:
the first --individuals whose two haplotypes each have the chromosome as a
single sequence go into the collection, and the haplotype 1 of the next
--queries such individuals become queries. A chromosome is found through
the UCSC chromAlias file of the haplotype's GenBank assembly, and only its
bases are downloaded, using the bgzip index (.fai and .gzi) of the
assembly on S3. Each sequence is upper-cased, characters other than
A/C/G/T become N, and it is reverse-complemented if most of its sampled
32-mers occur in the reverse complement of the reference.

Writes OUTDIR/collection.txt (the reference first, then the collection,
one sequence per line), OUTDIR/queries.txt, OUTDIR/ref_length.txt (the
length of the reference plus one, for build_index --reference-length) and
OUTDIR/manifest.tsv.
"""
import argparse
import csv
import io
import os
import random
import struct
import sys
import time
import urllib.request
import zlib

INDEX_URL = ("https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/"
             "main/data_tables/assemblies_release2_v1.0.index.csv")
UCSC_HUBS = "https://hgdownload.soe.ucsc.edu/hubs"
S3_HTTPS = "https://s3-us-west-2.amazonaws.com/"
COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")


def _acgt_table():
    """Upper-cases a, c, g and t and maps every other byte to N."""
    table = bytearray(b"N" * 256)
    for c in b"ACGT":
        table[c] = c
        table[c + 32] = c  # lower case
    return bytes(table)


ACGT = _acgt_table()


def https(url):
    return S3_HTTPS + url[len("s3://"):] if url.startswith("s3://") else url


def fetch(url, start=None, end=None, tries=4):
    """The bytes of url (a URL or a local path), or of its range [start, end)."""
    if "://" not in url:
        with open(url, "rb") as f:
            if start is None:
                return f.read()
            f.seek(start)
            return f.read(None if end is None else end - start)
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
        except Exception as e:  # network errors: retry
            if attempt + 1 == tries:
                raise
            print(f"  retrying {url}: {e}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))


def read_fai(url):
    fai = {}
    for line in fetch(url).decode().splitlines():
        name, length, offset, linebases, linewidth = line.split("\t")[:5]
        fai[name] = (int(length), int(offset), int(linebases), int(linewidth))
    return fai


def read_gzi(url):
    data = fetch(url)
    (n,) = struct.unpack_from("<Q", data, 0)
    pairs = [(0, 0)] + [struct.unpack_from("<QQ", data, 8 + 16 * k) for k in range(n)]
    return pairs  # (compressed offset, uncompressed offset) of each block start


def extract(fasta_url, fai, gzi, name):
    """The sequence `name` from a bgzipped FASTA, downloading only its blocks."""
    length, offset, linebases, linewidth = fai[name]
    u0 = offset
    u1 = offset + (length // linebases) * linewidth + length % linebases
    first = max(k for k in range(len(gzi)) if gzi[k][1] <= u0)
    later = [k for k in range(len(gzi)) if gzi[k][1] >= u1]
    c0, base = gzi[first]
    c1 = gzi[later[0]][0] if later else None
    data = fetch(fasta_url, c0, c1)
    out = io.BytesIO()
    while data:
        d = zlib.decompressobj(31)
        out.write(d.decompress(data))
        data = d.unused_data
    text = out.getvalue()[u0 - base:u1 - base]
    seq = text.replace(b"\n", b"").replace(b"\r", b"")
    if len(seq) != length:
        raise RuntimeError(f"{name}: extracted {len(seq)} bases, expected {length}")
    return seq.translate(ACGT)


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


class Orienter:
    """Decides whether a sequence matches the reference or its reverse
    complement, from 32-mers of the reference sampled every 64 positions."""

    def __init__(self, ref, k=32, step=64):
        self.k, self.step = k, step
        self.kmers = {ref[i:i + k] for i in range(0, len(ref) - k, step) if b"N" not in ref[i:i + k]}

    def hits(self, seq, rng, samples=2000):
        n = 0
        for _ in range(samples):
            i = rng.randrange(0, max(1, len(seq) - self.k - self.step))
            if any(seq[i + d:i + d + self.k] in self.kmers for d in range(self.step)):
                n += 1
        return n

    def orient(self, seq, seed):
        fwd = self.hits(seq, random.Random(seed))
        rc_seq = reverse_complement(seq)
        rev = self.hits(rc_seq, random.Random(seed))
        return (seq, "+", fwd, rev) if fwd >= rev else (rc_seq, "-", fwd, rev)


def chrom_alias_name(genbank_accession, chrom):
    """The PanSN name of `chrom` in a GenBank assembly, from its UCSC chromAlias
    file, or None if the chromosome is missing or in several pieces."""
    acc = genbank_accession  # e.g. GCA_041900255.1
    digits = acc[4:].split(".")[0]
    url = f"{UCSC_HUBS}/GCA/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/{acc}/{acc}.chromAlias.txt"
    try:
        text = fetch(url).decode()
    except Exception as e:
        print(f"  no chromAlias for {acc}: {e}", file=sys.stderr)
        return None
    header = None
    names = []
    for line in text.splitlines():
        if line.startswith("#"):
            header = line.lstrip("# ").split("\t")
            continue
        fields = dict(zip(header, line.split("\t")))
        if fields.get("ucsc") == chrom:
            names.append(fields.get("hprcV2"))
    return names[0] if len(names) == 1 else None


class Assembly:
    def __init__(self, row):
        self.row = row
        self.fasta = https(row["assembly"])
        self._fai = self._gzi = None

    def fai(self):
        if self._fai is None:
            self._fai = read_fai(https(self.row["assembly_fai"]))
        return self._fai

    def get(self, name):
        if self._gzi is None:
            self._gzi = read_gzi(https(self.row["assembly_gzi"]))
        return extract(self.fasta, self.fai(), self._gzi, name)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--chrom", default="chr19")
    ap.add_argument("--individuals", type=int, default=16)
    ap.add_argument("--queries", type=int, default=2)
    ap.add_argument("--index", default=INDEX_URL)
    ap.add_argument("--min-fraction", type=float, default=0.9,
                    help="minimum length of a haplotype's chromosome, as a fraction of the reference's")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = list(csv.DictReader(io.StringIO(fetch(args.index).decode())))
    chm13 = next(r for r in rows if r["sample_id"] == "CHM13")
    ref = Assembly(chm13).get(f"CHM13#0#{args.chrom}")
    print(f"reference CHM13 {args.chrom}: {len(ref)} bases", file=sys.stderr)
    orienter = Orienter(ref)

    by_sample = {}
    for r in rows:
        if r["haplotype"] in ("1", "2") and r["genbank_accession"].startswith("GCA_"):
            by_sample.setdefault(r["sample_id"], {})[r["haplotype"]] = r

    manifest = []
    collection, queries = [], []
    for sample in sorted(by_sample):
        if len(collection) >= 2 * args.individuals and len(queries) >= args.queries:
            break
        role = "collection" if len(collection) < 2 * args.individuals else "query"
        haps = ["1", "2"] if role == "collection" else ["1"]
        if any(h not in by_sample[sample] for h in haps):
            continue
        got = []
        for h in haps:
            row = by_sample[sample][h]
            name = chrom_alias_name(row["genbank_accession"], args.chrom)
            if name is None:
                break
            asm = Assembly(row)
            if name not in asm.fai() or asm.fai()[name][0] < args.min_fraction * len(ref):
                break
            seq, strand, fwd, rev = orienter.orient(asm.get(name), seed=len(manifest) + 1)
            if max(fwd, rev) < 0.5 * 2000:
                print(f"  {sample} hap{h}: too few shared 32-mers ({fwd}, {rev})", file=sys.stderr)
                break
            got.append((h, name, seq, strand, fwd, rev))
        if len(got) < len(haps):
            print(f"{sample}: skipped", file=sys.stderr)
            continue
        for h, name, seq, strand, fwd, rev in got:
            (collection if role == "collection" else queries).append(seq)
            manifest.append((role, sample, h, name, len(seq), seq.count(b"N"), strand, fwd, rev))
            print(f"{sample} hap{h} -> {role}: {name}, {len(seq)} bases, strand {strand}", file=sys.stderr)

    if len(collection) < 2 * args.individuals or len(queries) < args.queries:
        sys.exit(f"found only {len(collection)} collection haplotypes and {len(queries)} queries")
    with open(os.path.join(args.outdir, "collection.txt"), "wb") as f:
        for seq in [ref] + collection:
            f.write(seq + b"\n")
    with open(os.path.join(args.outdir, "queries.txt"), "wb") as f:
        for seq in queries:
            f.write(seq + b"\n")
    with open(os.path.join(args.outdir, "ref_length.txt"), "w") as f:
        f.write(f"{len(ref) + 1}\n")
    with open(os.path.join(args.outdir, "manifest.tsv"), "w") as f:
        f.write("role\tsample\thaplotype\tsequence\tlength\tN\tstrand\tkmer_hits_fwd\tkmer_hits_rc\n")
        f.write("reference\tCHM13\t0\tCHM13#0#%s\t%d\t%d\t+\t\t\n" % (args.chrom, len(ref), ref.count(b"N")))
        for m in manifest:
            f.write("\t".join(map(str, m)) + "\n")


if __name__ == "__main__":
    main()
