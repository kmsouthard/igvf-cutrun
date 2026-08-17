#!/usr/bin/env python3
"""Check every processing claim made about the Hs27 CUT&RUN signal files.

The claims are the ones submitted to IGVF -- assembly, deduplication, bin size,
normalization -- and each is re-derived here from the deposited bigWigs rather
than taken from the methods text or the filenames. Anything that cannot be
derived from a bigWig is reported as UNVERIFIABLE rather than passed.

Dependency-free: standard library only. Exits non-zero if any check fails.

    python3 verify.py <directory of the eight .bw files>
    python3 verify.py --fetch <work dir>     # pull them from Zenodo first

Zenodo record 15215216, Hs27_CnR_Density_bigwigs_hg38_10mNorm.zip.
"""
import hashlib
import os
import struct
import sys
import urllib.request
import zipfile
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bw_quantum as bq
from bw_summary import read_bw

ZENODO = ("https://zenodo.org/api/records/15215216/files/"
          "Hs27_CnR_Density_bigwigs_hg38_10mNorm.zip/content")

# --- the claims, as constants -------------------------------------------------
ASSEMBLY_CHR1 = 248_956_422          # GRCh38; mm10 chr1 is 195,471,971
CONTIGS = 194                        # no alt loci, no chrEBV
GENOME_BP = 3_099_750_718
BIN_BP = 25
NORM_TARGET = 10_000_000
ENTITY_BP_RANGE = (400, 600)         # what the observed means imply
BLOCKS = 300                         # bigWig data blocks to sample per file

# accession, md5 and raw depth are from the IGVF portal, not from these files
LIBRARIES = {
    "RA480-1": dict(target="H3K4me3",  rep=1, signal_file="IGVFFI4865LLEJ",
                    md5="b2c1016d03a3e5c409a9702f87b9624e", read_pairs=6315051),
    "RA480-2": dict(target="H3K4me3",  rep=2, signal_file="IGVFFI9030WSCX",
                    md5="250c127cbea75d021dfcba843a4418e7", read_pairs=5769753),
    "RA480-3": dict(target="H3K27me3", rep=1, signal_file="IGVFFI3631QNSO",
                    md5="c2f8174c2f786b2966b43a839a0c1da7", read_pairs=7515513),
    "RA480-4": dict(target="H3K27me3", rep=2, signal_file="IGVFFI1865NTKM",
                    md5="870496dce611679d5b09e89c77c53ad5", read_pairs=7685945),
    "RA480-5": dict(target="H3K27ac",  rep=1, signal_file="IGVFFI1472HREP",
                    md5="feec3f221a418eea1035f46cb819bcb2", read_pairs=7673653),
    "RA480-6": dict(target="H3K27ac",  rep=2, signal_file="IGVFFI7555NGXG",
                    md5="400f5a76aa3d48d140fb020897c99df1", read_pairs=6771430),
    "RA480-7": dict(target="IgG",      rep=1, signal_file="IGVFFI6525PXBV",
                    md5="3739b832c723d848bfb7d75206af78af", read_pairs=5064103),
    "RA480-8": dict(target="IgG",      rep=2, signal_file="IGVFFI3258CYNL",
                    md5="5c5da05e06cc81d4a3e90e3c05433f02", read_pairs=6276001),
}


class Report:
    def __init__(self):
        self.passed = self.failed = self.unverifiable = 0

    def ok(self, what, detail=""):
        self.passed += 1
        print(f"  PASS  {what}" + (f"  [{detail}]" if detail else ""))

    def bad(self, what, detail=""):
        self.failed += 1
        print(f"  FAIL  {what}" + (f"  [{detail}]" if detail else ""))

    def na(self, what, why):
        self.unverifiable += 1
        print(f"  n/a   {what}  [{why}]")

    def check(self, cond, what, detail=""):
        (self.ok if cond else self.bad)(what, detail)
        return cond


def fetch(workdir):
    os.makedirs(workdir, exist_ok=True)
    zp = os.path.join(workdir, "Hs27_CnR.zip")
    if not os.path.exists(zp):
        print(f"fetching {ZENODO}")
        urllib.request.urlretrieve(ZENODO, zp)
    with zipfile.ZipFile(zp) as z:
        for n in z.namelist():
            if n.endswith(".bw") and not os.path.basename(n).startswith("._"):
                z.extract(n, workdir)
    for root, _d, fs in os.walk(workdir):
        if any(f.endswith(".bw") for f in fs):
            return root
    raise SystemExit("no bigWigs in the archive")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def scan(path):
    """Bin grid, value quantum and totalSummary, in one pass over sampled blocks."""
    d = read_bw(path)
    sizes = dict(d["chroms"])
    by_id = {i: n for i, (n, _) in enumerate(sorted(d["chroms"]))}

    fh, endian, meta = bq._open(path)
    with fh:
        items = bq.leaves(fh, endian, meta["index"])
        step = max(1, len(items) // BLOCKS)
        values, starts_off, ends_off, ends_at_contig = [], 0, 0, 0
        for off, size in items[::step][:BLOCKS]:
            fh.seek(off)
            raw = fh.read(size)
            if meta["uncompress"]:
                raw = zlib.decompress(raw)
            cid, _a, _b, _s, _p, btype, _r, n = struct.unpack_from(
                endian + "IIIIIBBH", raw, 0)
            if btype != 1:                       # not bedGraph; grid is undefined
                continue
            clen = sizes.get(by_id.get(cid, ""), 0)
            body = raw[24:]
            for i in range(n):
                a, b, v = struct.unpack_from(endian + "IIf", body, i * 12)
                values.append(v)
                if a % BIN_BP:
                    starts_off += 1
                if b % BIN_BP:
                    ends_off += 1
                    if b == clen:
                        ends_at_contig += 1

    s = d["summary"]
    scale = bq.quantum(values)
    pos = [v for v in values if v > 0]
    exact = sum(1 for v in pos
                if abs(v / scale - round(v / scale)) < 2e-4 * max(round(v / scale), 1))
    return dict(sizes=sizes, mean=s["total"] / s["valid"], valid=s["valid"],
                scale=scale, n=len(values), n_pos=len(pos), exact=exact,
                starts_off=starts_off, ends_off=ends_off,
                ends_at_contig=ends_at_contig)


def main():
    argv = sys.argv[1:]
    if argv and argv[0] == "--fetch":
        bwdir = fetch(argv[1] if len(argv) > 1 else "zenodo")
    elif argv:
        bwdir = argv[0]
    else:
        raise SystemExit(__doc__)

    files = {}
    for f in sorted(os.listdir(bwdir)):
        if not f.endswith(".bw") or f.startswith("._"):
            continue
        for lib in LIBRARIES:
            if f"_{lib}_" in f:
                files[lib] = os.path.join(bwdir, f)

    r = Report()
    print(f"\n{len(files)} of {len(LIBRARIES)} bigWigs found in {bwdir}\n")
    r.check(len(files) == len(LIBRARIES), "all eight signal files present")

    rows = []
    for lib in sorted(LIBRARIES):
        meta = LIBRARIES[lib]
        print(f"\n{lib}  {meta['target']} rep{meta['rep']}  ->  {meta['signal_file']}")
        if lib not in files:
            r.bad("file present")
            continue
        path = files[lib]

        r.check(md5(path) == meta["md5"],
                "byte-identical to the file on the IGVF portal", meta["md5"])

        m = scan(path)
        r.check(m["sizes"].get("chr1") == ASSEMBLY_CHR1,
                "assembly is GRCh38", f"chr1 = {m['sizes'].get('chr1'):,}")
        r.check(len(m["sizes"]) == CONTIGS and sum(m["sizes"].values()) == GENOME_BP,
                "contig set is GRCh38 primary, no alts, no chrEBV",
                f"{len(m['sizes'])} contigs, {sum(m['sizes'].values()):,} bp")
        r.check(all(c.startswith("chr") for c in m["sizes"]),
                "all contigs are UCSC-named")
        r.check(m["valid"] == GENOME_BP,
                "coverage spans the whole reference, nothing blacklisted out")

        r.check(m["starts_off"] == 0,
                f"every interval starts on the {BIN_BP} bp grid",
                f"{m['n']:,} intervals")
        r.check(m["ends_off"] == m["ends_at_contig"],
                "ragged intervals occur only at contig ends",
                f"{m['ends_off']} ragged, {m['ends_at_contig']} at a contig end")

        r.check(m["exact"] == m["n_pos"],
                "all values are integer multiples of one scale factor",
                f"{m['n_pos']:,} positive values")

        depth = NORM_TARGET / m["scale"]
        keep = 100.0 * depth / meta["read_pairs"]
        r.check(0 < keep < 100,
                f"a {NORM_TARGET // 10**6}M target implies a possible dedup rate",
                f"scale {m['scale']:.2f}, {depth:,.0f} of {meta['read_pairs']:,} = {keep:.1f}%")
        keep2 = 100.0 * (2 * NORM_TARGET / m["scale"]) / meta["read_pairs"]

        entity = m["mean"] * GENOME_BP / NORM_TARGET - BIN_BP
        r.check(ENTITY_BP_RANGE[0] <= entity <= ENTITY_BP_RANGE[1],
                "mean implies fragment-sized, not read-sized, counting units",
                f"{entity:.0f} bp; 101 bp reads would give mean "
                f"{NORM_TARGET * (101 + BIN_BP) / GENOME_BP:.2f}, observed {m['mean']:.2f}")

        r.na("duplicates were removed",
             "a bigWig cannot show it; RmDup in the filename, confirmed at the bench")
        r.na("no MAPQ filter was applied",
             "not derivable; decides read-depth signal vs signal of unique reads")

        rows.append((lib, m["scale"], depth, keep, keep2, entity))

    print("\n" + "=" * 78)
    print("Normalization target: the dedup rate a target implies must be under 100%")
    print(f"{'library':10}{'scale':>7}{'  10M implies':>16}{'keep':>8}{'  20M keep':>11}")
    for lib, scale, depth, keep, keep2, _e in rows:
        flag = "  <- impossible" if keep2 >= 100 else ""
        print(f"{lib:10}{scale:>7.2f}{depth:>16,.0f}{keep:>7.1f}%{keep2:>10.1f}%{flag}")
    if rows:
        r.check(any(k2 >= 100 for *_x, k2, _e in rows),
                "\na 20M target is excluded outright by at least one library")

    print("\n" + "=" * 78)
    print(f"{r.passed} passed, {r.failed} failed, {r.unverifiable} not verifiable "
          f"from the files")
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
