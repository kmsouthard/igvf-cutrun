#!/usr/bin/env python3
"""Read the chromosome list out of a bigWig/bigBed header without pyBigWig.

Parses the B+ tree at chromosomeTreeOffset per the Kent bigWig spec.
"""
import struct
import sys

BIGWIG_MAGIC = 0x888FFC26
BIGBED_MAGIC = 0x8789F2EB
BTREE_MAGIC = 0x78CA8C91


def read_chroms(path):
    with open(path, "rb") as fh:
        raw = fh.read(4)
        if len(raw) < 4:
            raise ValueError("file too short")
        for endian in ("<", ">"):
            magic = struct.unpack(endian + "I", raw)[0]
            if magic in (BIGWIG_MAGIC, BIGBED_MAGIC):
                break
        else:
            raise ValueError("not a bigWig/bigBed (bad magic)")
        kind = "bigWig" if magic == BIGWIG_MAGIC else "bigBed"

        fh.seek(0)
        hdr = fh.read(64)
        (_magic, version, zoom_levels, chrom_tree_offset, full_data_offset,
         full_index_offset, field_count, defined_field_count, autosql_offset,
         total_summary_offset, uncompress_buf_size, _reserved) = struct.unpack(
            endian + "IHHQQQHHQQIQ", hdr)

        fh.seek(chrom_tree_offset)
        bmagic, block_size, key_size, val_size, item_count, _res = struct.unpack(
            endian + "IIIIQQ", fh.read(32))
        if bmagic != BTREE_MAGIC:
            raise ValueError("bad B+ tree magic at chromosomeTreeOffset")

        chroms = []

        def walk(offset):
            fh.seek(offset)
            is_leaf, _pad, count = struct.unpack(endian + "BBH", fh.read(4))
            if is_leaf:
                for _ in range(count):
                    key = fh.read(key_size).rstrip(b"\x00").decode(
                        "utf-8", "replace")
                    chrom_id, chrom_size = struct.unpack(
                        endian + "II", fh.read(val_size))
                    chroms.append((chrom_id, key, chrom_size))
            else:
                children = []
                for _ in range(count):
                    fh.read(key_size)
                    children.append(struct.unpack(endian + "Q", fh.read(8))[0])
                for child in children:
                    walk(child)

        walk(chrom_tree_offset + 32)
        chroms.sort()
        return kind, version, item_count, chroms


if __name__ == "__main__":
    for path in sys.argv[1:]:
        try:
            kind, version, n, chroms = read_chroms(path)
        except Exception as exc:  # noqa: BLE001 - report and continue
            print(f"{path}\n  ERROR: {exc}\n")
            continue
        print(f"{path}")
        print(f"  format={kind} version={version} n_chroms={n}")
        names = [c[1] for c in chroms]
        print(f"  first 30: {names[:30]}")
        prefixed = [x for x in names if x.startswith("chr")]
        print(f"  chr-prefixed: {len(prefixed)}/{len(names)}")
        odd = [x for x in names if not x.startswith("chr")]
        if odd:
            print(f"  NON chr-prefixed ({len(odd)}): {odd[:40]}")
        print()
