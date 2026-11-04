#!/usr/bin/env python3
"""The chromosome 19 data of "More RLBWT Tricks" v2 (tool/CHR19.md, "Data"), in the format
experiments/chr19_paper.sh reads: OUTDIR/collection.txt (one sequence per line, CHM13 first, then both
haplotypes of each individual) and OUTDIR/queries.txt (the two query haplotypes).

    python3 hprc_chr19.py OUTDIR [--manifest FILE] [--individuals N] [--plan]

The sequences come from HPRC release 2, as in chr19/twolevel/download_hprc_chr19.py, whose download code
this reuses (only the bases of chromosome 19 are fetched, block by block from the bgzipped assemblies on S3;
about 60 MB per sequence).  Every character other than A, C, G and T becomes N, and each sequence is
checked against its length and number of Ns.

--manifest FILE: a TSV with columns role (reference, collection or query), sample, haplotype, sequence
name in the assembly, length, number of N; the collection's individuals in order.  Without it, the 16
individuals and 2 queries of download_hprc_chr19.py.  --individuals N: the first N individuals of the
manifest (2N + 1 sequences with CHM13).  --plan: list the sequences and stop.
"""
import argparse, importlib.util, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(HERE, "..", "..", "chr19", "twolevel", "download_hprc_chr19.py")


def load_original():
    sys.dont_write_bytecode = True                    # (no __pycache__ in chr19/, which stays unchanged)
    spec = importlib.util.spec_from_file_location("download_hprc_chr19", ORIG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_manifest(path):
    rows = []
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        role, sample, hap, name, length, ns = line.rstrip("\n").split("\t")[:6]
        rows.append((role, sample, hap, name, int(length), int(ns)))
    return rows


def select(rows, individuals):
    ref = [r for r in rows if r[0] == "reference"]
    queries = [r for r in rows if r[0] == "query"]
    order = []
    for r in rows:
        if r[0] == "collection" and r[1] not in order:
            order.append(r[1])
    if individuals is not None:
        if individuals > len(order):
            sys.exit("the manifest has only %d individuals" % len(order))
        order = order[:individuals]
    coll = [r for s in order for r in rows if r[0] == "collection" and r[1] == s]
    for s in order:
        if sum(1 for r in coll if r[1] == s) != 2:
            sys.exit("individual %s does not have two haplotypes in the manifest" % s)
    if len(ref) != 1 or len(queries) != 2:
        sys.exit("the manifest needs one reference and two queries")
    return ref + coll, queries


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--manifest")
    ap.add_argument("--individuals", type=int)
    ap.add_argument("--plan", action="store_true")
    a = ap.parse_args()
    orig = load_original()
    rows = read_manifest(a.manifest) if a.manifest else list(orig.DATASET)
    coll, queries = select(rows, a.individuals)
    print("%d sequences in the collection (%d individuals), %d queries" % (len(coll), (len(coll) - 1) // 2, len(queries)), file=sys.stderr)
    if a.plan:
        for r in coll + queries:
            print("\t".join(map(str, r)))
        return
    os.makedirs(a.outdir, exist_ok=True)
    index = {}
    import csv, io
    for r in csv.DictReader(io.StringIO(orig.fetch(orig.INDEX_URL).decode())):
        index[(r["sample_id"], r["haplotype"])] = r
    parts = os.path.join(a.outdir, "parts"); os.makedirs(parts, exist_ok=True)
    for k, (role, sample, hap, name, length, ns) in enumerate(coll + queries):
        path = os.path.join(parts, name.replace("#", "_") + ".seq")
        if os.path.exists(path) and os.path.getsize(path) == length:
            continue
        print("[%d/%d] %s (%s)" % (k + 1, len(coll) + len(queries), name, role), file=sys.stderr)
        seq = orig.extract(index[(sample, hap)], name).translate(orig.ACGT)
        if seq.count(b"N") != ns:
            sys.exit("%s: %d Ns, expected %d" % (name, seq.count(b"N"), ns))
        with open(path + ".tmp", "wb") as f:
            f.write(seq)
        os.replace(path + ".tmp", path)
    for stem, entries in (("collection", coll), ("queries", queries)):
        with open(os.path.join(a.outdir, stem + ".txt"), "wb") as out:
            for r in entries:
                with open(os.path.join(parts, r[3].replace("#", "_") + ".seq"), "rb") as f:
                    out.write(f.read() + b"\n")
    print("wrote %s/collection.txt and queries.txt (the per-sequence files are kept in %s for reruns)" % (a.outdir, parts), file=sys.stderr)


if __name__ == "__main__":
    main()
