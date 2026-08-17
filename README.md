# igvf-cutrun

Processing pipeline for the Hs27 fibroblast CUT&RUN dataset submitted to the IGVF data
portal, from *Comprehensive transcription factor perturbations recapitulate fibroblast
transcriptional states* (Nat Genet 2025,
[doi:10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)).

Eight libraries — two biological replicates each of H3K4me3, H3K27me3, H3K27ac and an
IgG control — from parental Hs27 fibroblasts, paired-end 101 bp on a NovaSeq 6000. The
submitted output is one depth-normalized bigWig per library.

| Target | Analysis set | Measurement set | Libraries |
|---|---|---|---|
| H3K4me3 | [IGVFDS5465EDQZ](https://data.igvf.org/analysis-sets/IGVFDS5465EDQZ/) | [IGVFDS6141QFVR](https://data.igvf.org/measurement-sets/IGVFDS6141QFVR/) | `RA480-1`, `RA480-2` |
| H3K27me3 | [IGVFDS8522VRST](https://data.igvf.org/analysis-sets/IGVFDS8522VRST/) | [IGVFDS8014KPTJ](https://data.igvf.org/measurement-sets/IGVFDS8014KPTJ/) | `RA480-3`, `RA480-4` |
| H3K27ac | [IGVFDS6938GALL](https://data.igvf.org/analysis-sets/IGVFDS6938GALL/) | [IGVFDS7154HHWN](https://data.igvf.org/measurement-sets/IGVFDS7154HHWN/) | `RA480-5`, `RA480-6` |
| IgG (control) | [IGVFDS0089EIMW](https://data.igvf.org/analysis-sets/IGVFDS0089EIMW/) | [IGVFDS3000EQWJ](https://data.igvf.org/measurement-sets/IGVFDS3000EQWJ/) | `RA480-7`, `RA480-8` |

`pipeline.sh` is a reconstruction: the libraries were processed by MSKCC's Epigenomics
Research Innovation Lab and their script did not survive, so the commands below are
rebuilt from the published methods and from the deposited files, which `verify/` checks
they reproduce.

## The pipeline

Reference throughout is GRCh38 primary assembly in UCSC contig naming — 194 sequences,
3,099,750,718 bp, no alt loci and no `chrEBV`
([`GCA_000001405.15_GRCh38_no_alt_analysis_set`](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/),
IGVF `IGVFFI0653VCGH`).

```
./pipeline.sh RA480-1 reads/RA480-1_R1.fastq.gz reads/RA480-1_R2.fastq.gz
```

### 1 · Trim — Trim Galore! 0.4.5

Adapter and quality trimming at Q15, calling cutadapt 1.15 and FastQC 0.11.5.

```bash
trim_galore --paired --quality 15 --fastqc --output_dir trimmed "$R1" "$R2"
```

### 2 · Align — bowtie2 2.3.4.1

```bash
bowtie2 --threads 8 -x "$REF" -1 "$T1" -2 "$T2" \
  | samtools sort -o "bam/$LIB.sorted.bam"
```

### 3 · Deduplicate — Picard 2.16.0

```bash
java -jar picard.jar MarkDuplicates \
     INPUT="bam/$LIB.sorted.bam" OUTPUT="bam/$LIB.sorted.RmDup.bam" \
     METRICS_FILE="bam/$LIB.dup_metrics.txt" REMOVE_DUPLICATES=true
```

Retains 32–50% of read pairs across the eight libraries.

### 4 · Call peaks — MACS2 2.1.2

IgG as the control, then blacklist filtering. Peak files are *not* part of the IGVF
submission; the step is here because it is part of the published analysis.

```bash
macs2 callpeak -t "bam/$LIB.sorted.RmDup.bam" -c bam/RA480-7.sorted.RmDup.bam \
               -f BAMPE -g hs -p 0.001 -n "$LIB" --outdir peaks
bedtools intersect -v -a "peaks/${LIB}_peaks.narrowPeak" -b "$BLACKLIST" \
       > "peaks/$LIB.peaks.bed"
```

### 5 · Normalized signal — BEDTools 2.29.2 → bedGraphToBigWig

**This step produces the submitted files.** Fragments are counted in fixed 25 bp
windows and scaled to a depth of 10 million, then written as bigWig. The blacklist is
not applied here — coverage spans the whole reference.

```bash
N=$(( $(samtools view -c -f 0x2 -F 0x100 "bam/$LIB.sorted.RmDup.bam") / 2 ))
SCALE=$(python3 -c "print(f'{10000000/$N:.2f}')")

bedtools makewindows -g "$CHROM_SIZES" -w 25 > signal/windows_25bp.bed

bedtools bamtobed -bedpe -i "bam/$LIB.sorted.RmDup.bam" \
  | awk 'BEGIN{OFS="\t"} $1==$4 {print $1, ($2<$5?$2:$5), ($3>$6?$3:$6)}' \
  | sort -k1,1 -k2,2n > "signal/$LIB.fragments.bed"

bedtools intersect -c -a signal/windows_25bp.bed -b "signal/$LIB.fragments.bed" \
  | awk -v s="$SCALE" 'BEGIN{OFS="\t"} {print $1, $2, $3, $4*s}' \
  | merge_equal_runs > "signal/$LIB.bedGraph"

bedGraphToBigWig "signal/$LIB.bedGraph" "$CHROM_SIZES" \
                 "signal/${LIB}_hg38_sorted.RmDup.10mNorm.bw"
```

Scale factors run 2.58 to 5.57, one per library, each rounded to two decimals —
`results/bigwig_forensics.tsv` has them all.

### 6 · Average replicates — deepTools 3.5.4

```bash
bigwigAverage -b signal/RA480-1_*.bw signal/RA480-2_*.bw -o signal/H3K4me3_mean.bw
```

Averaged tracks are *not* part of the IGVF submission; the eight files there are all
per-replicate.

## Verifying

`verify/verify.py` fetches the deposited bigWigs from Zenodo, checks their md5sums
against the IGVF portal records, then re-derives every processing claim from the bytes.
Standard library only; exits non-zero on failure.

```
python3 verify/verify.py --fetch work/
python3 verify/verify.py path/to/Density_bigwigs_hg38_10mNorm/    # local copy
```

Current result: **82 passed, 0 failed, 16 not verifiable from the files.**

| Step | Claim | How it is established |
|---|---|---|
| — | identity | md5 matches the IGVF record for that accession |
| 2 | GRCh38 | `chr1` is 248,956,422 — mm10's is 195,471,971 |
| 2 | no-alt reference | 194 UCSC-named contigs, 3,099,750,718 bp |
| 3 | duplicates removed | **not verifiable** — reported `n/a`, not passed |
| 3 | no MAPQ filter | **not verifiable** — reported `n/a`, not passed |
| 4 | blacklist not applied to coverage | intervals span the whole reference, nothing zeroed |
| 5 | 25 bp windows | every interval start divides by 25; ragged intervals only at contig ends |
| 5 | one scale factor per file | 100% of positive values are exact integer multiples of the smallest |
| 5 | target was 10 million | `10e6 / scale` against portal `read_count` gives 32–50% duplicate retention; 20 million would force `RA480-4` above 100% and is excluded |
| 5 | fragments, not reads | in 25 bp bins an entity of length `L` occupies `(L+25)/25` bins, so mean is `1e7(L+25)/genome`. 101 bp reads predict 0.41; observed is 1.52–1.74, giving `L` of 445–515 bp |

The last two are why this is code rather than a note: both are non-obvious, and the
fragments-not-reads result is not what the published wording says.

`PIPELINE.md` marks every parameter above as `[methods]`, `[measured]` or `[open]`, and
lists the nine still outstanding.

## Layout

```
pipeline.sh    the pipeline, end to end
PIPELINE.md    every parameter, with its source and what is still open
verify/        verify.py, plus three standalone bigWig readers it is built from
results/       measured per-file output, and the saved verify run
metadata/      the portal records these measurements imply, as submission sheets
```

`verify/bw_summary.py`, `bw_chroms.py` and `bw_quantum.py` work on any bigWig.
`bw_quantum.py` recovers a normalization scale factor from the quantization of the
stored values, which is how the 10 million target was established independently of the
filename.

## Related

* Paper: [doi:10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)
  · preprint [PMC11312553](https://pmc.ncbi.nlm.nih.gov/articles/PMC11312553/)
* Data: Zenodo [10.5281/zenodo.15215216](https://doi.org/10.5281/zenodo.15215216),
  `Hs27_CnR_Density_bigwigs_hg38_10mNorm.zip` · raw reads SRA
  [PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1108254)
* IGVF portal: [Tom Norman lab, MSKCC](https://data.igvf.org/labs/tom-norman/) —
  analysis sets are visible once released
* Companion repositories: `igvf-atac`, `igvf-rnaseq`, and
  [kmsouthard/igvf-crispra](https://github.com/kmsouthard/igvf-crispra) for the
  Perturb-seq element-level analysis
