# Hs27 CUT&RUN — processing specification

What produced the eight signal files submitted to IGVF. Sources are marked throughout:
**[methods]** is quoted or paraphrased from the published paper, **[measured]** is
re-derived from the deposited bigWigs by `verify/verify.py`, and **[open]** is
neither — a parameter the files constrain but the methods do not state.

The published source is the Methods section *"Epigenetic data processing and
preparation of model features"*. Read it as covering ATAC and CUT&RUN together; it
does.

## Inputs

Eight libraries, `RA480-1` to `RA480-8`: two biological replicates each of H3K4me3,
H3K27me3, H3K27ac and IgG, from parental Hs27 fibroblasts. Paired-end 101 bp on an
Illumina NovaSeq 6000, all eight on one run and one lane, 5.1–7.7 M read pairs each.
Accessions and per-library read counts are in `verify/verify.py`.

## Steps

### 1. Read trimming

**[methods]** Raw reads were trimmed and filtered for quality (Q>15) and adapter
content using Trim Galore! 0.4.5, running cutadapt 1.15 and FastQC 0.11.5.

**[measured]** The FASTQs deposited on the portal are uniformly 101 bp, so this step
operates downstream of what was submitted — which is correct, the raw reads are the
submitted artifact.

### 2. Alignment

**[methods]** bowtie2 2.3.4.1, to human assembly hg38.

**[measured]** GRCh38 confirmed: `chr1` is 248,956,422. The contig set is 194
sequences in UCSC naming — 25 primary, 42 unlocalized `*_random`, 127 unplaced
`chrUn_*` — totalling 3,099,750,718 bp. No alt loci and no `chrEBV`.

**[open]** Which hg38 FASTA. The measured contig set is the no-alt primary assembly,
not the full 455-sequence UCSC `hg38.fa`. It matches
`GCA_000001405.15_GRCh38_no_alt_analysis_set` to within the single `chrEBV` decoy
contig, which is why IGVF `IGVFFI0653VCGH` is used as `reference_files` on the
submitted files.

**[open]** `samtools` is absent from the methods but `sorted` is in every output
filename, so BAM sorting sits here and needs a version.

### 3. Duplicate removal

**[methods]** Alignments were deduplicated using MarkDuplicates in Picard Tools 2.16.0
(RRID:SCR_006525).

**[not verifiable]** A bigWig cannot show whether deduplication ran. The claim rests on
`RmDup` in every filename plus confirmation from the lab, and the harness reports it as
`n/a` rather than passing it.

**[open]** Whether a mapping-quality or unique-mapper filter was also applied. This is
the only open item that changes a submitted field: without one the files are IGVF
`content_type: read-depth signal`, which is what was submitted; with one they would be
`signal of unique reads`.

### 4. Peak calling

**[methods]** Enriched regions were called with MACS2 2.1.2 at p 0.001 and filtered
against a blacklist.

**[measured, consistent]** The blacklist was applied to peaks only, not to coverage:
the bigWigs span the whole reference with nothing zeroed.

Peak files exist but are **not** part of the IGVF submission; the four analysis sets
contain only the eight bigWigs. This step is recorded for completeness.

### 5. Signal generation and normalization

The step that produced the submitted files.

**[methods]** "The BEDTools suite v2.29.2 was used to create bigwig files normalized to
10 million mapped reads."

**[measured]** Binned at exactly 25 bp. Every interval start in every file divides by
25 — 745,492 of 745,492 in the largest sample — and the only ragged intervals end
exactly at a contig length, i.e. the truncated final window of each contig.

**[measured]** One global scale factor per file: 100% of positive values are exact
integer multiples of that file's smallest positive value. Factors run 2.58 to 5.57,
each rounded to two decimals, which indicates a factor computed outside the counting
tool and passed in.

**[measured]** The 10 million target is confirmed independently of the methods text.
`10,000,000 / scale` gives the deduplicated depth implied; against the portal's
`read_count` for each library that must be under 100%. A 10M target gives 32–50%
retention across all eight; a 20M target forces `RA480-4` to 100.9% and is excluded
outright, putting a hard ceiling at 19.8M.

**[open]** The exact invocation. `bedtools genomecov` is base-pair resolution and
cannot produce a 25 bp grid; `bedtools makewindows -w 25` followed by `coverage`
produces exactly it, terminal truncation included. The attribution to BEDTools is
consistent, the description is incomplete, and 25 bp is a real parameter of the
submitted files.

**[open]** BEDTools writes bedGraph, not bigWig, so a `bedGraphToBigWig` step exists
and needs a version.

**[open]** What the 10 million counted. Read literally as individual reads it does not
reproduce the files. Counting in 25 bp bins, an entity of length `L` occupies
`(L + 25) / 25` bins, so a 10M-normalized track has genome-wide mean
`1e7 x (L + 25) / genome`, independent of library depth. For 101 bp reads that predicts
a mean near **0.41**; the observed means are **1.52 to 1.74**, requiring `L` between
**445 and 515 bp**. So the denominator counted fragments — read pairs — or reads were
extended to fragment length. This halves or doubles every implied depth. It changes no
submitted field, and the phrase "10 million mapped reads" is kept verbatim from the
methods in the file descriptions.

### 6. Replicate averaging

**[methods]** Normalized bigWigs were averaged across replicates with deepTools 3.5.4
`bigwigAverage`.

The averaged tracks are **not** part of the IGVF submission; the eight submitted files
are all per-replicate. Recorded for completeness.

## Software

Versions as used. Portal status checked 2026-08-17 against 251 registered Software
objects: no tool is registered at the version this pipeline used, so every
SoftwareVersion is new even where the Software object exists.

| Step | Tool | Version | IGVF Software object |
|---|---|---|---|
| 1 | Trim Galore! | 0.4.5 | none — create |
| 1 | cutadapt | 1.15 | `cutadapt`, has 3.4.0 only |
| 1 | FastQC | 0.11.5 | `fastqc`, has 0.11.9 only |
| 2 | bowtie2 | 2.3.4.1 | `bowtie2`, has 2.3.4.3 |
| 2 | samtools | **[open]** | `samtools` |
| 3 | Picard Tools | 2.16.0 | none — create |
| 4 | MACS2 | 2.1.2 | `macs`, has 2.1.0 / 2.2.4 / 2.2.7.1 |
| 5 | BEDTools | 2.29.2 | `bedtools`, has 2.29.0 / 2.30.0 |
| 5 | bedGraphToBigWig | **[open]** | `bedgraphtobigwig` |
| 6 | deepTools | 3.5.4 | none — create |

Only steps 1, 2, 3 and 5 are needed for what is currently submitted.

## What carries over to `igvf-atac`

**ATAC-seq** is covered by the same methods paragraph — steps 1 through 5 apply
unchanged, and the open questions have the same answers. One thing does not carry over:
a **second, independent** ATAC processing run exists, nf-core/atacseq 2.1.2, whose full
`software_versions.yml` survives. It is a different pipeline with a different
normalization (CPM rather than 10M) and a different contig set (146 or 166 rather than
194). Which pipeline produced the *submitted* ATAC tracks must be settled before either
set of versions is used, and reading the contig count off those files answers it.

**Bulk RNA-seq** shares nothing here. It is nf-core/rnaseq with STAR and salmon, and
`igvf-rnaseq` needs its own specification from scratch.
