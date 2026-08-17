#!/usr/bin/env python3
"""Identify a bigWig's genome assembly and summarise its signal, header-only.

Dependency-free. Two things a submission needs and a filename cannot prove:

  assembly  -- inferred from the primary chromosome LENGTHS, which are unique
               per build (chr1 is 248,956,422 in GRCh38 and 249,250,621 in hg19).
  signal    -- from the bigWig totalSummary block: coverage breadth, mean, and
               dynamic range. Enough to tell a focal mark from a broad one, a
               normalised track from raw counts, and a real IP from a control.

Usage:
    python3 bw_summary.py file1.bw file2.bw ...
"""
import struct
import sys
import os

BIGWIG_MAGIC = 0x888FFC26
BTREE_MAGIC = 0x78CA8C91

# chr1 length -> assembly. Unique per build, so one contig identifies it.
CHR1_LENGTH = {
    248956422: "GRCh38 / hg38",
    249250621: "GRCh37 / hg19",
    247249719: "NCBI36 / hg18",
    195471971: "GRCm38 / mm10",
    197195432: "NCBI37 / mm9",
    195154279: "GRCm39 / mm39",
}


def read_bw(path):
    with open(path, "rb") as fh:
        raw = fh.read(4)
        for endian in ("<", ">"):
            if struct.unpack(endian + "I", raw)[0] == BIGWIG_MAGIC:
                break
        else:
            raise ValueError("not a bigWig")

        fh.seek(0)
        hdr = struct.unpack(endian + "IHHQQQHHQQIQ", fh.read(64))
        chrom_tree_offset, total_summary_offset = hdr[3], hdr[9]

        # chromosome B+ tree
        fh.seek(chrom_tree_offset)
        bmagic, _bs, key_size, val_size, _n, _r = struct.unpack(
            endian + "IIIIQQ", fh.read(32))
        if bmagic != BTREE_MAGIC:
            raise ValueError("bad chrom tree")
        chroms = []

        def walk(off):
            fh.seek(off)
            is_leaf, _pad, count = struct.unpack(endian + "BBH", fh.read(4))
            if is_leaf:
                for _ in range(count):
                    key = fh.read(key_size).rstrip(b"\x00").decode("utf-8", "replace")
                    cid, size = struct.unpack(endian + "II", fh.read(val_size))
                    chroms.append((key, size))
            else:
                kids = []
                for _ in range(count):
                    fh.read(key_size)
                    kids.append(struct.unpack(endian + "Q", fh.read(8))[0])
                for k in kids:
                    walk(k)

        walk(chrom_tree_offset + 32)

        summary = None
        if total_summary_offset:
            fh.seek(total_summary_offset)
            valid, mn, mx, total, sq = struct.unpack(endian + "Qdddd", fh.read(40))
            summary = dict(valid=valid, min=mn, max=mx, total=total, sumsq=sq)
        return dict(chroms=chroms, summary=summary)


def main():
    for path in sys.argv[1:]:
        print(os.path.basename(path))
        try:
            d = read_bw(path)
        except Exception as exc:  # noqa: BLE001
            print(f"  ERROR: {exc}\n")
            continue

        sizes = dict(d["chroms"])
        asm = "unknown"
        for name in ("chr1", "1"):
            if name in sizes:
                asm = CHR1_LENGTH.get(sizes[name], f"unrecognised (chr1={sizes[name]:,})")
                break
        genome_bp = sum(sizes.values())
        print(f"  assembly     : {asm}")
        print(f"  contigs      : {len(sizes)}  ({sum(1 for c in sizes if c.startswith('chr'))} chr-prefixed)")
        print(f"  chrY present : {'yes' if any(c in sizes for c in ('chrY','Y')) else 'no'}")

        s = d["summary"]
        if s and s["valid"]:
            mean = s["total"] / s["valid"]
            var = max(s["sumsq"] / s["valid"] - mean * mean, 0.0)
            sd = var ** 0.5
            breadth = 100.0 * s["valid"] / genome_bp if genome_bp else float("nan")
            print(f"  covered bp   : {s['valid']:,}  ({breadth:.1f}% of contigs in file)")
            print(f"  mean signal  : {mean:.4f}   sd {sd:.4f}")
            print(f"  min / max    : {s['min']:.4f} / {s['max']:.2f}")
            print(f"  max/mean     : {s['max']/mean:,.0f}x   <- focal marks run high, controls low")
            intish = abs(mean - round(mean)) < 1e-9 and abs(s["max"] - round(s["max"])) < 1e-9
            print(f"  values       : {'integer-like (raw counts?)' if intish else 'non-integer (scaled/normalised)'}")
        print()


if __name__ == "__main__":
    main()
