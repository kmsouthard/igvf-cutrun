#!/usr/bin/env bash
#
# Hs27 CUT&RUN: FASTQ to normalized bigWig.
#
# Steps 2-4 are the command lines recorded in the delivered files' own headers.
# Step 5 has no surviving invocation and is reconstructed. See PIPELINE.md.
#
#     ./pipeline.sh RA480-1 reads/RA480-1_R1.fastq.gz reads/RA480-1_R2.fastq.gz
#
set -euo pipefail

LIB=${1:?library, e.g. RA480-1}
R1=${2:?R1 fastq}
R2=${3:?R2 fastq}

REF=ref/GCA_000001405.15_GRCh38_no_alt_analysis_set     # bowtie2 index prefix
CHROM_SIZES=ref/hg38.no_alt.chrom.sizes                 # 194 contigs, 3,099,750,718 bp
# Peak-calling control. The delivered peaks used a shared human CD34+ ChIP input
# [recorded]; IgG is the matched control this assay actually provides. Override to
# swap, or set to "none" to call without one. The output name follows the control,
# so files stay self-describing.
CONTROL=${CONTROL:-ref/humanCD34_input.hg38.sorted.RmDup.bam}
BIN=25                                                  # [reconstructed]
MINMAPQ=10                                              # [reconstructed]
TARGET=10000000                                         # 10 million, per the methods
THREADS=${THREADS:-8}

mkdir -p trimmed bam peaks signal

# --- 1. Trim and QC ---------------------------------------- Trim Galore! 0.4.5
#     invokes cutadapt 1.15 and FastQC 0.11.5. Only the output filenames are
#     recorded; the parameters below are from the methods.
trim_galore --paired --quality 15 --fastqc --cores 4 \
            --output_dir trimmed "$R1" "$R2"

T1=trimmed/$(basename "${R1%.fastq.gz}")_val_1.fq.gz
T2=trimmed/$(basename "${R2%.fastq.gz}")_val_2.fq.gz

# --- 2. Align ------------------------------------------------------ bowtie2 2.3.4.1
#     [recorded] --local -p 6. Local mode, not end-to-end.
bowtie2 --local -p 6 -x "$REF" -1 "$T1" -2 "$T2" \
  | samtools sort -@ "$THREADS" -o "bam/$LIB.hg38.sorted.bam"
samtools index "bam/$LIB.hg38.sorted.bam"

# --- 3. Remove duplicates ------------------------------- Picard 2.16.0, MarkDuplicates
#     [recorded] REMOVE_DUPLICATES=true ASSUME_SORTED=true. No MAPQ filter here.
java -jar picard.jar MarkDuplicates \
     INPUT="bam/$LIB.hg38.sorted.bam" \
     OUTPUT="bam/$LIB.hg38.sorted.RmDup.bam" \
     METRICS_FILE="bam/$LIB.dup_metrics.txt" \
     REMOVE_DUPLICATES=true \
     ASSUME_SORTED=true
samtools index "bam/$LIB.hg38.sorted.RmDup.bam"

# --- 4. Call peaks --------------------------------------------------- MACS2 2.1.2
#     [recorded] -f BAM (not BAMPE), --nomodel, p 0.001. As delivered, every
#     library including IgG was called against the CD34+ input.
TREAT="bam/$LIB.hg38.sorted.RmDup.bam"

case ${CONTROL:-none} in
  none|"")
    CTL_ARGS=(); CTL_TAG=noControl ;;
  *humanCD34_input*)
    CTL_ARGS=(-c "$CONTROL"); CTL_TAG=vsInput ;;              # as delivered
  *)
    CTL_ARGS=(-c "$CONTROL")
    CTL_TAG="vs$(basename "$CONTROL" | cut -d. -f1)" ;;
esac

# -ef is same-device-and-inode, and is false if either file is absent
if [[ ${CONTROL:-none} != none && "$CONTROL" -ef "$TREAT" ]]; then
  echo "$LIB: control is the treatment; skipping peak call" >&2
else
  macs2 callpeak -t "$TREAT" ${CTL_ARGS[@]+"${CTL_ARGS[@]}"} \
                 -f BAM -g hs --nomodel --pvalue 0.001 \
                 -n "MACS2_${CTL_TAG}_$LIB.hg38.sorted.RmDup" --outdir peaks
fi

# --- 5. Normalized signal --------------------------------------- [reconstructed]
#     No invocation survives. The fitted parameters -- 25 bp bins, fragments,
#     MAPQ >= 10, duplicates already gone -- are what bamCoverage produces, and a
#     sibling project's log records exactly these. Scale is 1e7 / fragments,
#     rounded to two decimals.
#
#     Fragments = pairs with both mates mapped at MAPQ >= 10, counted once. This
#     definition reproduces the deposited scale factors to within 0.1-2.7%;
#     proper-pair-only and mapped/2 do not. verify/fragcount.py computes it.
N=$(samtools view -c -f 0x40 -F 0x90c -q "$MINMAPQ" "bam/$LIB.hg38.sorted.RmDup.bam")
SCALE=$(python3 -c "print(f'{$TARGET/$N:.2f}')")
echo "$LIB: $N fragments, scale $SCALE"

bamCoverage --bam "bam/$LIB.hg38.sorted.RmDup.bam" \
            --outFileName "signal/${LIB}_hg38_sorted.RmDup.10mNorm.bw" \
            --outFileFormat bigwig \
            --binSize "$BIN" \
            --scaleFactor "$SCALE" \
            --minMappingQuality "$MINMAPQ" \
            --ignoreDuplicates \
            --extendReads \
            --maxPairedFragmentLength 1000 \
            --numberOfProcessors "$THREADS"

# --- 6. Average replicates ------------------------------- deepTools 3.5.4
#     Once per target, over its two replicates.
#
#     bigwigAverage -b signal/RA480-1_*.bw signal/RA480-2_*.bw \
#                   -o signal/H3K4me3_mean.bw
