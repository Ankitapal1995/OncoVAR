#!/bin/bash

# Print usage
usage() {
    echo "Usage: $0 --ref <reference.gb> --R1 <reads_1.fastq.gz> --R2 <reads_2.fastq.gz> --prefix <prefix> --output <output_dir>"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --ref) REF="$2"; shift ;;
        --R1) READS1="$2"; shift ;;
        --R2) READS2="$2"; shift ;;
        --prefix) PREFIX="$2"; shift ;;
        --output) OUTPUT="$2"; shift ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# Check if mandatory arguments are provided
if [[ -z "$REF" || -z "$READS1" || -z "$READS2" || -z "$PREFIX" || -z "$OUTPUT" ]]; then
    echo "Error: Missing arguments."
    usage
fi
START_TIME=$(date +%s)

# Create output directory
mkdir -p "$OUTPUT"

# Run the pipeline
set -e


if [[ "$REF" == *.fna ]]; then
   BASE_NAME=$(basename "$REF" .gbk)
elif [[ "$REF" == *.gb ]]; then
    BASE_NAME=$(basename "$REF" .gb)
elif [["$REF" == *.fna ]]; then 
     BASE_NAME=$(basename "$REF" .fna)
else
    echo "Error: File must have .gb or .gbk extension."
    exit 1
fi


# Define the output file
REF_fasta="${BASE_NAME}"

echo "Running FastQC..."
fastqc "$READS1" "$READS2" -o "$OUTPUT"

echo "Running Trimmomatic..."
trimmomatic PE -phred33 \
    "$READS1" "$READS2" \
    "$OUTPUT/${PREFIX}_trimmed_R1.fastq.gz" "$OUTPUT/${PREFIX}_unpaired_R1.fastq.gz" \
    "$OUTPUT/${PREFIX}_trimmed_R2.fastq.gz" "$OUTPUT/${PREFIX}_unpaired_R2.fastq.gz" \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36 -threads 8

echo "Running BWA..."
bwa index "$REF_fasta"
bwa mem -t 8 -T 0 -M -R "@RG\tID:${PREFIX}\tLB:${PREFIX}\tPL:ILLUMINA\tPM:HISEQ\tSM:${PREFIX}" \
    "$REF_fasta" "$OUTPUT/${PREFIX}_trimmed_R1.fastq.gz" "$OUTPUT/${PREFIX}_trimmed_R2.fastq.gz" > "$OUTPUT/${PREFIX}_aligned_reads.sam"

echo "Converting SAM to BAM and sorting..."
samtools view -S -b "$OUTPUT/${PREFIX}_aligned_reads.sam" > "$OUTPUT/${PREFIX}_aligned.bam"
samtools sort "$OUTPUT/${PREFIX}_aligned.bam" -o "$OUTPUT/${PREFIX}_sorted.bam"
samtools sort -n "$OUTPUT/${PREFIX}_sorted.bam" -o "$OUTPUT/${PREFIX}_query_sorted.bam"

echo "Running Samtools Fixmate..."
samtools fixmate -m "$OUTPUT/${PREFIX}_query_sorted.bam" "$OUTPUT/${PREFIX}_fixed.bam"


samtools sort "$OUTPUT/${PREFIX}_fixed.bam" -o "$OUTPUT/${PREFIX}_sorted_fixed.bam"

echo "Marking duplicates..."
samtools markdup -r "$OUTPUT/${PREFIX}_sorted_fixed.bam" "$OUTPUT/${PREFIX}_deduplicated.bam"

echo "Indexing deduplicated bam file"
samtools index "$OUTPUT/${PREFIX}_deduplicated.bam"

echo "Pipeline completed successfully. Results are in the $OUTPUT directory."
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Total time taken: $DURATION seconds"
