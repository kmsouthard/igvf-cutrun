# Parameters

`[methods]` is from the published Methods section *"Epigenetic data processing and
preparation of model features"*. `[recorded]` is read out of the delivered files' own
headers — bowtie2 and Picard `@PG` lines in the BAMs, the MACS2 command line in each
`_peaks.xls`. `[measured]` is re-derived from the deposited bigWigs by `verify/verify.py`.
`[open]` is none of these.

Source of the `[recorded]` values: `results/RA480_pipeline_headers.txt`, extracted from the
ETIL `rerunAug2026` delivery. Identical across all eight libraries.

## 1 · Trim

* `[methods]` Trim Galore! 0.4.5, cutadapt 1.15, FastQC 0.11.5, quality Q > 15.
* `[recorded]` Trim Galore ran: bowtie2's inputs are `*_val_1.fq.gz` / `*_val_2.fq.gz`.
* `[open]` Trim Galore's own parameters. Only its output filenames survive.

## 2 · Align

* `[methods]` bowtie2 2.3.4.1, human assembly hg38.
* `[recorded]` `bowtie2-align-s --local -p 6 -x <bowtie2-index>/hg38`. **Local mode**,
  not end-to-end — the methods do not state this.
* `[recorded]` 194 contigs, 3,099,750,718 bp in the BAM `@SQ` block; `@HD SO:coordinate`.
* `[measured]` The same 194-contig set, independently, from the bigWigs.
* `[open]` Which hg38 FASTA built the bowtie2 index. The contig set is the no-alt primary
  assembly, not the 455-sequence UCSC `hg38.fa`.
* `[open]` samtools version. No samtools `@PG` line survives, though the BAM is sorted.

## 3 · Deduplicate

* `[methods]` Picard 2.16.0 `MarkDuplicates`.
* `[recorded]` `MarkDuplicates VN:2.16.0-1-g763d98e-SNAPSHOT`, `REMOVE_DUPLICATES=true`,
  `ASSUME_SORTED=true`, `DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES`,
  `OPTICAL_DUPLICATE_PIXEL_DISTANCE=100`. Confirms the published 2.16.0 for this dataset.
* `[recorded]` No MAPQ filter at this step — the BAM retains MAPQ 0 alignments.

## 4 · Call peaks

* `[methods]` MACS2 2.1.2, p 0.001, blacklist filtered.
* `[recorded]` `macs2 callpeak -t <lib>.bam -c <control> -f BAM -g hs --nomodel
  --pvalue 0.001`. Effective genome 2.70e9, tag size 96 bp, `d = 200`,
  **Paired-End mode off** — `-f BAM`, not `-f BAMPE`, so mates enter as independent tags.
* `[recorded]` **The control is not IgG.** All eight libraries were called against
  `humanCD34_input.hg38.sorted.RmDup.bam`, a generic human CD34+
  ChIP input from the core's shared resources (55.4M tags). The two IgG libraries were
  themselves called as treatment against it.
* `[recorded]` No blacklist step appears in the MACS2 invocation; the delivered
  `_peaks.narrowPeak` files are raw MACS2 output. Any blacklist filtering the methods
  describe happened downstream, or not at all.
* `pipeline.sh` takes the control as `$CONTROL`, defaulting to the CD34+ input so the
  default run reproduces the delivered filenames. Point it at an IgG library, or set it to
  `none`, to call peaks the other ways; the `MACS2_<tag>_` prefix follows the control, and a
  library is skipped rather than called against itself.
* `[open]` Which control the November 2023 peak calls used. Peak counts differ from the
  2023 `Total_number_peaks.xlsx` by +9.1% from slightly *fewer* usable fragments, which
  implicates the control rather than the alignment — see
  `results/rerun2026_vs_2023_consistency.tsv`. A contemporaneous log for a sibling project
  names a different control,
  `Input.hg38.sorted.RmDup.bam`.

## 5 · Normalized signal

The one step with no surviving invocation. The parameters below are reconstructed.

* `[methods]` BEDTools 2.29.2, normalized to 10 million mapped reads.
* `[measured]` Binned at 25 bp: every interval start divides by 25, and the only ragged
  intervals end at a contig length.
* `[measured]` One global scale factor per file; 100% of positive values are exact integer
  multiples of it. Factors 2.58–5.57, two decimals — see `results/bigwig_forensics.tsv`.
* `[measured]` Target is 10 million.
* `[measured]` The counting unit is the **fragment**. In 25 bp bins an entity of length `L`
  occupies `(L + 25) / 25` bins, so the genome-wide mean is `1e7 (L + 25) / genome`. For
  101 bp reads that predicts 0.41; observed means are 1.52–1.74, giving `L` = 445–515 bp.
* `[measured]` **A MAPQ >= 10 filter was applied.** Counting each BAM as pairs with both
  mates mapped at MAPQ >= 10, once per pair, reproduces `1e7 / scale` to within
  −0.1% .. −2.7% (six of eight within 1%). Rival definitions do not fit: proper-pair-only
  −3.5% .. −12.9%, mapped/2 +6.8% .. +11.8%, pairs-plus-singletons +17% .. +37%. Check it
  with `verify/fragcount.py`. This is a fit, not a logged parameter.
* `[open]` The exact invocation. The fitted parameters — 25 bp bins, fragments, MAPQ >= 10,
  duplicates already removed — are what `deeptools bamCoverage` produces, and a sibling
  project's log records exactly `binLength 25`, `minMappingQuality 10`,
  `ignoreDuplicates True`, `maxPairedFragmentLength 1000`. BEDTools cannot write bigWig, so
  the published attribution to BEDTools 2.29.2 is doubtful for this step.
* `[open]` `bedGraphToBigWig` version, if it was used at all.

## 6 · Average replicates

* `[methods]` deepTools 3.5.4 `bigwigAverage`.

## Software

| Step | Tool | Version | Source |
|---|---|---|---|
| 1 | Trim Galore! · cutadapt · FastQC | 0.4.5 · 1.15 · 0.11.5 | `[methods]` |
| 2 | bowtie2 | 2.3.4.1 | `[recorded]` |
| 2 | samtools | `[open]` | |
| 3 | Picard Tools | 2.16.0 | `[recorded]` |
| 4 | MACS2 | 2.1.2 | `[recorded]` |
| 5 | deepTools *or* BEDTools | `[open]` | |
| 5 | bedGraphToBigWig | `[open]` | |
| 6 | deepTools | 3.5.4 | `[methods]` |

Steps 1, 2, 3 and 5 produced the bigWigs deposited at IGVF. Alignments and peaks were
delivered by the core in August 2026 and are not yet deposited.
