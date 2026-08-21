#!/usr/bin/env python3
"""Count fragments under deepTools bamCoverage filters, straight from BAM records.
Filters mirror the KS467 process log: minMappingQuality 10, ignoreDuplicates,
maxPairedFragmentLength 1000, properly-paired, one count per pair."""
import gzip, struct, sys, os

def count(path, mapq_min=10, max_frag=1000):
    f = gzip.open(path, 'rb')
    assert f.read(4) == b'BAM\x01'
    l_text = struct.unpack('<i', f.read(4))[0]; f.read(l_text)
    n_ref = struct.unpack('<i', f.read(4))[0]
    for _ in range(n_ref):
        ln = struct.unpack('<i', f.read(4))[0]; f.read(ln); f.read(4)

    tot = 0            # all records
    mapped = 0         # not unmapped, not secondary/supplementary
    frag_all = 0       # proper pair, first-in-pair, primary  (no MAPQ, no len cut)
    frag_mapq = 0      # + MAPQ>=10
    frag_full = 0      # + |tlen| <= max_frag
    unp = struct.unpack_from
    buf = b''; pos = 0
    CH = 1 << 22
    while True:
        if len(buf) - pos < 4:
            buf = buf[pos:] + f.read(CH); pos = 0
            if len(buf) < 4: break
        bs = unp('<i', buf, pos)[0]; pos += 4
        if len(buf) - pos < bs:
            buf = buf[pos:] + f.read(max(CH, bs)); pos = 0
            if len(buf) < bs: break
        mapq = buf[pos + 9]
        flag = unp('<H', buf, pos + 14)[0]
        tlen = unp('<i', buf, pos + 28)[0]
        pos += bs
        tot += 1
        if flag & 0x904:                      # unmapped | secondary | supplementary
            continue
        mapped += 1
        if (flag & 0x2) and (flag & 0x40):    # proper pair, first in pair
            frag_all += 1
            if mapq >= mapq_min:
                frag_mapq += 1
                if abs(tlen) <= max_frag:
                    frag_full += 1
    f.close()
    return tot, mapped, frag_all, frag_mapq, frag_full

if __name__ == '__main__':
    for p in sys.argv[1:]:
        t, m, fa, fq, ff = count(p)
        print(f'{os.path.basename(p)}\t{t}\t{m}\t{fa}\t{fq}\t{ff}', flush=True)
