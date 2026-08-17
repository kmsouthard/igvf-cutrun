#!/usr/bin/env bash
#
# Hs27 CUT&RUN: FASTQ to normalized bigWig.
#
# Versions and parameters are from the published methods; lines marked
# [reconstructed] are recovered from the deposited files. See PIPELINE.md.
#
#     ./pipeline.sh RA480-1 reads/RA480-1_R1.fastq.gz reads/RA480-1_R2.fastq.gz
#
set -euo pipefail

LIB=${1:?library, e.g. RA480-1}
R1=${2:?R1 fastq}
R2=${3:?R2 fastq}

REF=ref/GCA_000001405.15_GRCh38_no_alt_analysis_set     # bowtie2 index prefix
CHROM_SIZES=ref/hg38.no_alt.chrom.sizes                 # 194 contigs, 3,099,750,718 bp
BLACKLIST=ref/hg38-blacklist.bed
BIN=25                                                  # [reconstructed]
TARGET=10000000                                         # 10 million, per the methods
THREADS=${THREADS:-8}

mkdir -p trimmed bam peaks signal

# --- 1. Trim and QC ---------------------------------------- Trim Galore! 0.4.5
#     invokes cutadapt 1.15 and FastQC 0.11.5
trim_galore --paired --quality 15 --fastqc --cores 4 \
            --output_dir trimmed "$R1" "$R2"

T1=trimmed/$(basename "${R1%.fastq.gz}")_val_1.fq.gz
T2=trimmed/$(basename "${R2%.fastq.gz}")_val_2.fq.gz

# --- 2. Align ------------------------------------------------------ bowtie2 2.3.4.1
#     GRCh38 primary assembly in UCSC naming, no alt loci, no chrEBV
bowtie2 --threads "$THREADS" -x "$REF" -1 "$T1" -2 "$T2" \
  | samtools sort -@ "$THREADS" -o "bam/$LIB.sorted.bam"
samtools index "bam/$LIB.sorted.bam"

# --- 3. Remove duplicates ------------------------------- Picard 2.16.0, MarkDuplicates
java -jar picard.jar MarkDuplicates \
     INPUT="bam/$LIB.sorted.bam" \
     OUTPUT="bam/$LIB.sorted.RmDup.bam" \
     METRICS_FILE="bam/$LIB.dup_metrics.txt" \
     REMOVE_DUPLICATES=true
samtools index "bam/$LIB.sorted.RmDup.bam"

# --- 4. Call peaks --------------------------------------------------- MACS2 2.1.2
#     IgG (RA480-7, RA480-8) is the control for the six histone libraries.
if [[ $LIB != RA480-7 && $LIB != RA480-8 ]]; then
  macs2 callpeak -t "bam/$LIB.sorted.RmDup.bam" -c bam/RA480-7.sorted.RmDup.bam \
                 -f BAMPE -g hs -p 0.001 -n "$LIB" --outdir peaks
  bedtools intersect -v -a "peaks/${LIB}_peaks.narrowPeak" -b "$BLACKLIST" \
         > "peaks/$LIB.peaks.bed"
fi

# --- 5. Normalized signal ------------------------- BEDTools 2.29.2, bedGraphToBigWig
#     Fragments counted in fixed 25 bp windows, scaled to 10 million. The
#     blacklist is not applied here.

# fragment count: properly paired, primary alignments only, divided by two
N=$(( $(samtools view -c -f 0x2 -F 0x100 "bam/$LIB.sorted.RmDup.bam") / 2 ))
SCALE=$(python3 -c "print(f'{$TARGET/$N:.2f}')")        # [reconstructed] 2 dp rounding
echo "$LIB: $N fragments, scale $SCALE"

bedtools makewindows -g "$CHROM_SIZES" -w "$BIN" > "signal/windows_${BIN}bp.bed"

bedtools bamtobed -bedpe -i "bam/$LIB.sorted.RmDup.bam" \
  | awk 'BEGIN{OFS="\t"} $1==$4 {print $1, ($2<$5?$2:$5), ($3>$6?$3:$6)}' \
  | sort -k1,1 -k2,2n > "signal/$LIB.fragments.bed"

# count per window, scale, then merge runs of equal value -- adjacent identical
# bins become one interval, as in the deposited files
bedtools intersect -c -a "signal/windows_${BIN}bp.bed" -b "signal/$LIB.fragments.bed" \
  | awk -v s="$SCALE" 'BEGIN{OFS="\t"} {print $1, $2, $3, $4*s}' \
  | awk 'BEGIN{OFS="\t"}
         NR>1 && $1==c && $2==e && $4==v {e=$3; next}
         NR>1 {print c, s, e, v}
         {c=$1; s=$2; e=$3; v=$4}
         END {print c, s, e, v}' \
  > "signal/$LIB.bedGraph"

bedGraphToBigWig "signal/$LIB.bedGraph" "$CHROM_SIZES" \
                 "signal/${LIB}_hg38_sorted.RmDup.10mNorm.bw"

# --- 6. Average replicates ------------------------------- deepTools 3.5.4
#     Once per target, over its two replicates.
#
#     bigwigAverage -b signal/RA480-1_*.bw signal/RA480-2_*.bw \
#                   -o signal/H3K4me3_mean.bw
