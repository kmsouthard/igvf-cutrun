# igvf-cutrun

CUT&RUN processing for Hs27 fibroblasts, from *Comprehensive transcription factor
perturbations recapitulate fibroblast transcriptional states* (Nat Genet 2025,
[doi:10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)).

Eight libraries — two replicates each of H3K4me3, H3K27me3, H3K27ac and IgG — to one
depth-normalized bigWig per library.

**`pipeline.sh` is a reconstruction.** These libraries were processed by MSKCC's
Epigenomics Research Innovation Lab and their script did not survive. The commands are
rebuilt from the published methods and from the deposited files; `PIPELINE.md` marks
every parameter as `[methods]`, `[measured]` or `[open]`, and `verify/` checks that the
reconstruction accounts for the files.

## Run

```bash
./pipeline.sh RA480-1 reads/RA480-1_R1.fastq.gz reads/RA480-1_R2.fastq.gz
```

Requires Trim Galore! 0.4.5, bowtie2 2.3.4.1, Picard 2.16.0, MACS2 2.1.2, BEDTools
2.29.2, bedGraphToBigWig and deepTools 3.5.4. Reference is GRCh38 primary assembly in
UCSC contig naming, no alt loci — 194 sequences, 3,099,750,718 bp
([`GCA_000001405.15_GRCh38_no_alt_analysis_set`](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/)).

Raw reads: SRA [PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1108254).

## Outputs

```
signal/<LIB>_hg38_sorted.RmDup.10mNorm.bw    coverage, 25 bp bins, scaled to 10M
peaks/<LIB>.peaks.bed                        MACS2 p 0.001, blacklist filtered
```

Per-library scale factors are in `results/bigwig_forensics.tsv`.

## Verify

```
python3 verify/verify.py --fetch work/                     # from Zenodo
python3 verify/verify.py path/to/Density_bigwigs_hg38_10mNorm/
```

Checks each deposited bigWig against its IGVF md5, then re-derives the assembly,
reference contig set, bin size, scale factor and normalization target from the bytes.
Deduplication and the absence of a MAPQ filter are reported `n/a` — a bigWig cannot
settle either. Currently 82 passed, 0 failed, 16 n/a.

`verify/bw_summary.py`, `bw_chroms.py` and `bw_quantum.py` work on any bigWig.

## Related

* Data: Zenodo [10.5281/zenodo.15215216](https://doi.org/10.5281/zenodo.15215216)
* [kmsouthard/igvf-atac](https://github.com/kmsouthard/igvf-atac) ·
  [kmsouthard/igvf-crispra](https://github.com/kmsouthard/igvf-crispra)
