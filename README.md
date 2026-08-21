# igvf-cutrun

CUT&RUN processing for Hs27 fibroblasts, from *Comprehensive transcription factor
perturbations recapitulate fibroblast transcriptional states* (Nat Genet 2025,
[doi:10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)).

Eight libraries — two replicates each of H3K4me3, H3K27me3, H3K27ac and IgG — to one
depth-normalized bigWig per library.

These libraries were processed by MSKCC's Epigenomics Research Innovation Lab. Align,
deduplicate and call-peaks are the command lines recorded in the delivered files' headers
(`results/RA480_pipeline_headers.txt`); signal generation has no surviving invocation and is
reconstructed. `PIPELINE.md` marks every parameter `[methods]`, `[recorded]`, `[measured]`
or `[open]` — including three that differ from the published methods.

## Run

```bash
./pipeline.sh RA480-1 reads/RA480-1_R1.fastq.gz reads/RA480-1_R2.fastq.gz
```

Requires Trim Galore! 0.4.5, bowtie2 2.3.4.1, samtools, Picard 2.16.0, MACS2 2.1.2 and
deepTools. Reference is GRCh38 primary assembly in UCSC contig naming, no alt loci — 194
sequences, 3,099,750,718 bp
([`GCA_000001405.15_GRCh38_no_alt_analysis_set`](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/)).

Raw reads: SRA [PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1108254).

## Outputs

```
bam/<LIB>.hg38.sorted.RmDup.bam                 duplicates removed, coordinate sorted
signal/<LIB>_hg38_sorted.RmDup.10mNorm.bw       coverage, 25 bp bins, scaled to 10M
peaks/MACS2_<tag>_<LIB>..._peaks.narrowPeak     MACS2 p 0.001, raw output; <tag> names the control
```

Scale factors: `results/bigwig_forensics.tsv`. August 2026 alignments and peaks against the
November 2023 delivery: `results/rerun2026_vs_2023_consistency.tsv`.

## Verify

```
python3 verify/verify.py --fetch work/                     # from Zenodo
python3 verify/verify.py path/to/Density_bigwigs_hg38_10mNorm/
python3 verify/fragcount.py bam/<LIB>.hg38.sorted.RmDup.bam
```

`verify.py` checks each deposited bigWig against its IGVF md5, then re-derives the assembly,
contig set, bin size, scale factor and normalization target from the bytes — 82 passed,
0 failed, 16 n/a. `fragcount.py` counts fragments under each candidate definition, to test
which reproduces a deposited scale factor. `bw_summary.py`, `bw_chroms.py` and
`bw_quantum.py` work on any bigWig.

## Data

Deposited at the IGVF Data Portal, one measurement set and one analysis set per target:

| Target | Measurement set | Analysis set |
|---|---|---|
| H3K4me3 | [IGVFDS6141QFVR](https://data.igvf.org/measurement-sets/IGVFDS6141QFVR/) | [IGVFDS5465EDQZ](https://data.igvf.org/analysis-sets/IGVFDS5465EDQZ/) |
| H3K27me3 | [IGVFDS8014KPTJ](https://data.igvf.org/measurement-sets/IGVFDS8014KPTJ/) | [IGVFDS8522VRST](https://data.igvf.org/analysis-sets/IGVFDS8522VRST/) |
| H3K27ac | [IGVFDS7154HHWN](https://data.igvf.org/measurement-sets/IGVFDS7154HHWN/) | [IGVFDS6938GALL](https://data.igvf.org/analysis-sets/IGVFDS6938GALL/) |
| IgG | [IGVFDS3000EQWJ](https://data.igvf.org/measurement-sets/IGVFDS3000EQWJ/) | [IGVFDS0089EIMW](https://data.igvf.org/analysis-sets/IGVFDS0089EIMW/) |

The bigWigs under `signal/` are the files in the analysis sets; the alignments and peaks are
not yet deposited. These records are `in progress`, so the links resolve only for signed-in
submitters until release.

## Related

* Data: Zenodo [10.5281/zenodo.15215216](https://doi.org/10.5281/zenodo.15215216)
* [kmsouthard/igvf-atac](https://github.com/kmsouthard/igvf-atac) ·
  [kmsouthard/igvf-crispra](https://github.com/kmsouthard/igvf-crispra)
