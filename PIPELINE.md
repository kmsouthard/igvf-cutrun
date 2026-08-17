# Parameters

`[methods]` is from the published Methods section *"Epigenetic data processing and
preparation of model features"*. `[measured]` is re-derived from the deposited bigWigs
by `verify/verify.py`. `[open]` is neither.

## 1 · Trim

* `[methods]` Trim Galore! 0.4.5, cutadapt 1.15, FastQC 0.11.5, quality Q > 15.

## 2 · Align

* `[methods]` bowtie2 2.3.4.1, human assembly hg38.
* `[measured]` GRCh38: `chr1` = 248,956,422. 194 contigs in UCSC naming — 25 primary,
  42 unlocalized, 127 unplaced — totalling 3,099,750,718 bp. No alt loci, no `chrEBV`.
* `[open]` Which hg38 FASTA. The measured set is the no-alt primary assembly, not the
  455-sequence UCSC `hg38.fa`.
* `[open]` samtools version. Absent from the methods; `sorted` is in every filename.

## 3 · Deduplicate

* `[methods]` Picard 2.16.0 `MarkDuplicates`.
* `[open]` Whether a MAPQ or unique-mapper filter also ran. Not derivable from a bigWig.

## 4 · Call peaks

* `[methods]` MACS2 2.1.2, p 0.001, blacklist filtered.
* `[measured]` The blacklist was applied to peaks only — coverage spans the whole
  reference with nothing zeroed.

## 5 · Normalized signal

* `[methods]` BEDTools 2.29.2, normalized to 10 million mapped reads.
* `[measured]` Binned at 25 bp: every interval start divides by 25, and the only ragged
  intervals end at a contig length.
* `[measured]` One global scale factor per file; 100% of positive values are exact
  integer multiples of it. Factors 2.58–5.57, two decimals — see
  `results/bigwig_forensics.tsv`.
* `[measured]` Target is 10 million. `1e7 / scale` against the portal's `read_count`
  gives 32–50% duplicate retention; 20 million forces `RA480-4` past 100% and is
  excluded.
* `[open]` The exact invocation. `bedtools genomecov` is base-pair resolution and cannot
  produce a 25 bp grid; `makewindows -w 25` plus `coverage` produces exactly it,
  terminal truncation included.
* `[open]` `bedGraphToBigWig` version. BEDTools writes bedGraph, not bigWig.
* `[open]` What the 10 million counted. In 25 bp bins an entity of length `L` occupies
  `(L + 25) / 25` bins, so the genome-wide mean is `1e7 (L + 25) / genome`. For 101 bp
  reads that predicts 0.41; observed means are 1.52–1.74, giving `L` = 445–515 bp. So
  fragments, or reads extended to fragment length.

## 6 · Average replicates

* `[methods]` deepTools 3.5.4 `bigwigAverage`.

## Software

| Step | Tool | Version |
|---|---|---|
| 1 | Trim Galore! · cutadapt · FastQC | 0.4.5 · 1.15 · 0.11.5 |
| 2 | bowtie2 | 2.3.4.1 |
| 2 | samtools | `[open]` |
| 3 | Picard Tools | 2.16.0 |
| 4 | MACS2 | 2.1.2 |
| 5 | BEDTools | 2.29.2 |
| 5 | bedGraphToBigWig | `[open]` |
| 6 | deepTools | 3.5.4 |

Steps 1, 2, 3 and 5 produced the files deposited at IGVF. Peaks and averaged tracks
were not deposited.
