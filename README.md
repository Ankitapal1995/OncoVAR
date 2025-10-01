# OncoVAR
OncoVAR is a custom pipeline for detecting cancer-specific somatic mutations from matched tumor-normal NGS data. It uses samtools and bcftools for variant calling, benchmarks against Mutect2 and VarScan2, and validates findings with a reference gene list to ensure accurate identification of tumor-exclusive mutations.


## TASK 2

### Requirement
Download `GATK` before proceeding.

### Steps

#### i) Quality Check, Alignment, and BAM File Preparation
Run the following command:
```bash
sh pipeline_for_bam_file.sh --ref <reference.gb> --R1 <reads_1.fastq.gz> --R2 <reads_2.fastq.gz> --prefix <prefix> --output <output_dir>
```

#### ii) Running Mutect2 and Varscan2
Use this command to perform mutation calling:
```bash
python mutect2_varscan2_mutation_calling.py reference.fasta tumor.bam normal.bam normal_sample_prefix
```

#### iii) Running In-House Pipeline
Execute the in-house somatic mutation calling pipeline with the following commands:
```bash
python pipeline_for_somatic_mut_calling.py reference.fasta normal.bam tumor.bam
python identify_cancer_somatic_mutation_cp.py tumor.vcf normal.vcf
```

#### iv) Background Estimation of GATK Mutect2-Generated Normal Sample VCF File
Use the following command for background estimation:
```bash
python background_estimation_new.py normal_sample_vcf_file total_bases output_prefix
```

To calculate the total bases, run:
```bash
grep -v '^>' GCF_000001405.40_GRCh38.p14_genomic.fna | tr -d '\n' | wc -c
```

