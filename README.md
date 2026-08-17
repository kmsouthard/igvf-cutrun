# igvf-cutrun

Processing provenance for the Hs27 fibroblast CUT&RUN dataset submitted to the IGVF
data portal, from *Comprehensive transcription factor perturbations recapitulate
fibroblast transcriptional states* (Nat Genet 2025,
`doi:10.1038/s41588-025-02284-1`).

Eight signal files across four analysis sets: H3K4me3 (`IGVFDS5465EDQZ`), H3K27ac
(`IGVFDS6938GALL`), H3K27me3 (`IGVFDS8522VRST`) and an IgG control
(`IGVFDS0089EIMW`), two biological replicates each.

Companion repositories cover the other assays from the same study: `igvf-atac`,
`igvf-rnaseq`, and [kmsouthard/igvf-crispra](https://github.com/kmsouthard/igvf-crispra)
for the Perturb-seq element-level analysis.

## What this repository is, and is not

**It is not the pipeline that produced the data.** These libraries were processed by
MSKCC's Epigenomics Research Innovation Lab using off-the-shelf tools. No lab-authored
script exists for that processing, and no run directory, workflow file or software
manifest from it survives.

**It is two things instead:**

1. **A specification** of what the pipeline did — tools, versions, parameters and
   reference — in `PIPELINE.md`, assembled from the published methods and completed
   where the methods are silent.
2. **A verification harness** — under `verify/` — that re-derives every processing
   claim made to IGVF directly from the deposited files, so each claim can be checked
   rather than taken on trust.

Anything the deposited files cannot settle is reported as `n/a`, never as a pass. Two
claims fall in that category and are labelled as such in the output: that duplicates
were removed, and that no mapping-quality filter was applied.

This distinction matters for anyone citing this repository as the `source_url` of an
IGVF AnalysisStepVersion. It records and checks the pipeline; it did not run it.

## Verifying

Standard library only, no dependencies. The harness fetches the deposited files from
Zenodo, checks their md5sums against the IGVF portal records, and then re-derives
everything from the bytes:

```
python3 verify/verify.py --fetch work/
```

or, against a local copy:

```
python3 verify/verify.py path/to/Density_bigwigs_hg38_10mNorm/
```

Exits non-zero if any check fails. As of the tag on this commit: **82 passed, 0 failed,
16 not verifiable from the files.**

### What it checks, per file

| Claim | How it is established |
|---|---|
| identity | md5 matches the IGVF portal record for that accession |
| GRCh38 | `chr1` is 248,956,422 — mm10's is 195,471,971 |
| reference contig set | 194 contigs, UCSC-named, 3,099,750,718 bp: primary assembly, no alt loci, no `chrEBV` |
| no blacklist filtering of coverage | intervals span the whole reference, nothing zeroed out |
| 25 bp bins | every interval start divides by 25; the only ragged intervals end exactly at a contig length |
| a single normalization factor | 100% of positive values are exact integer multiples of the file's smallest positive value |
| the target was 10 million | `target / scale` is the deduplicated depth implied; against portal `read_count` it must be under 100%. 10M gives 32–50%; 20M forces `RA480-4` to 100.9% and is excluded |
| counting units were fragments | in 25 bp bins an entity of length `L` occupies `(L+25)/25` bins, so a 10M track has mean `1e7(L+25)/genome`. 101 bp reads predict 0.41; observed means are 1.52–1.74, giving `L` of 445–515 bp |

The last two are the reason this harness exists rather than a prose note: both are
non-obvious, both are checkable, and one of them contradicts the plain reading of the
published methods.

## Layout

```
PIPELINE.md  what the pipeline did, step by step, with versions and open questions
verify/      the harness, plus three standalone bigWig readers it is built from
results/     measured per-file output of the harness
metadata/    the portal records these measurements imply, as submission sheets
```

`verify/bw_summary.py`, `verify/bw_chroms.py` and `verify/bw_quantum.py` are usable on
their own against any bigWig. `bw_quantum.py` recovers a normalization scale factor
from the quantization of the stored values, which is how the 10 million target was
established independently of the filename.

## Related

* Data: Zenodo [10.5281/zenodo.15215216](https://doi.org/10.5281/zenodo.15215216),
  `Hs27_CnR_Density_bigwigs_hg38_10mNorm.zip` · raw reads SRA `PRJNA1108254` ·
  IGVF portal, Tom Norman lab, MSKCC
* `igvf-atac` shares most of this pipeline — the published methods paragraph covers
  ATAC and CUT&RUN together. See the end of `PIPELINE.md` for what carries over and
  the one thing that does not.
