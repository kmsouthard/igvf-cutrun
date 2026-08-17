#!/usr/bin/env python3
"""Recover a bigWig's normalisation scale factor from the stored values.

A depth-normalised coverage track is raw integer coverage multiplied by a
single scale factor, so every value in the file is an integer multiple of that
factor. Recovering the factor recovers the normalisation:

    scale = 1 / (reads that produced one unit of coverage)
    implied library depth = target / scale

which is what a submission needs when the pipeline that made the file is not
otherwise documented. Dependency-free: parses the R-tree index and inflates the
data blocks directly.

Usage:
    python3 bw_quantum.py file.bw [file2.bw ...]
"""
import collections
import os
import struct
import sys
import zlib

BIGWIG_MAGIC = 0x888FFC26
RTREE_MAGIC = 0x2468ACE0
BLOCKS = 400          # leaf blocks to sample, spread across the file


def _open(path):
    fh = open(path, "rb")
    raw = fh.read(4)
    for endian in ("<", ">"):
        if struct.unpack(endian + "I", raw)[0] == BIGWIG_MAGIC:
            break
    else:
        fh.close()
        raise ValueError("not a bigWig")
    fh.seek(0)
    hdr = struct.unpack(endian + "IHHQQQHHQQIQ", fh.read(64))
    return fh, endian, dict(index=hdr[5], uncompress=hdr[10])


def leaves(fh, endian, index_offset):
    """Every leaf item in the R-tree, as (data_offset, data_size)."""
    fh.seek(index_offset)
    magic = struct.unpack(endian + "I", fh.read(4))[0]
    if magic != RTREE_MAGIC:
        raise ValueError("bad R-tree magic")
    fh.seek(index_offset + 48)                      # past the R-tree header

    out = []
    stack = [index_offset + 48]
    while stack:
        fh.seek(stack.pop())
        is_leaf, _pad, count = struct.unpack(endian + "BBH", fh.read(4))
        if is_leaf:
            buf = fh.read(count * 32)
            for i in range(count):
                off, size = struct.unpack_from(endian + "QQ", buf, i * 32 + 16)
                out.append((off, size))
        else:
            buf = fh.read(count * 24)
            for i in range(count):
                stack.append(struct.unpack_from(endian + "Q", buf, i * 24 + 16)[0])
    return out


def block_values(fh, endian, uncompress, off, size):
    fh.seek(off)
    raw = fh.read(size)
    if uncompress:
        raw = zlib.decompress(raw)
    _cid, _start, _end, _step, _span, btype, _r, n = struct.unpack_from(
        endian + "IIIIIBBH", raw, 0)
    body = raw[24:]
    if btype == 1:                                  # bedGraph
        return [struct.unpack_from(endian + "f", body, i * 12 + 8)[0]
                for i in range(n)]
    if btype == 2:                                  # varStep
        return [struct.unpack_from(endian + "f", body, i * 8 + 4)[0]
                for i in range(n)]
    if btype == 3:                                  # fixedStep
        return list(struct.unpack_from(endian + f"{n}f", body, 0))
    raise ValueError(f"unknown section type {btype}")


def quantum(values, tol=2e-4):
    """Largest q such that every value is a near-integer multiple of q.

    Starts from the smallest positive value -- coverage of 1 read -- and only
    shrinks if some value is not a multiple of it.
    """
    pos = sorted(v for v in values if v > 0)
    if not pos:
        return None
    q = pos[0]
    for v in pos:
        k = v / q
        if abs(k - round(k)) > tol * max(round(k), 1):
            # v is not a multiple of q; the true unit divides their difference
            q = min(q, abs(v - round(v / q) * q)) or q
    return q


def main():
    for path in sys.argv[1:]:
        print(os.path.basename(path))
        try:
            fh, endian, meta = _open(path)
        except Exception as exc:                     # noqa: BLE001
            print(f"  ERROR: {exc}\n")
            continue
        with fh:
            items = leaves(fh, endian, meta["index"])
            step = max(1, len(items) // BLOCKS)
            sample = items[::step][:BLOCKS]
            vals = []
            for off, size in sample:
                try:
                    vals.extend(block_values(fh, endian, meta["uncompress"], off, size))
                except Exception:                    # noqa: BLE001
                    continue

        pos = [v for v in vals if v > 0]
        if not pos:
            print("  no positive values in the sample\n")
            continue
        q = quantum(vals)
        mult = collections.Counter(round(v / q) for v in pos)
        exact = sum(1 for v in pos if abs(v / q - round(v / q)) < 2e-4 * max(round(v / q), 1))

        print(f"  blocks sampled  : {len(sample):,} of {len(items):,}")
        print(f"  values sampled  : {len(vals):,}  ({len(pos):,} positive)")
        print(f"  smallest > 0    : {min(pos):.8f}")
        print(f"  scale factor    : {q:.8f}   <- coverage of one read")
        print(f"  multiples exact : {100.0 * exact / len(pos):.2f}% of positive values")
        print(f"  implied depth   : {1e7 * q:,.0f} reads if the target was 10M")
        print(f"                    {1.0 / q:,.0f} reads if the target was 1x coverage")
        top = ", ".join(f"{k}x:{c:,}" for k, c in mult.most_common(6))
        print(f"  commonest units : {top}")
        print()


if __name__ == "__main__":
    main()
