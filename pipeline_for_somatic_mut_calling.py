import subprocess
import sys
def run_command(command):
    """Run a shell command and handle errors."""
    try:
        print(f"Running: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        exit(1)

def call_variants(input_bam, reference, output_vcf):
    """Call variants using bcftools."""
    command = (
        f"bcftools mpileup -f {reference} {input_bam} | "
        f"bcftools call -mv -Oz -o {output_vcf}"
    )
    run_command(command)

def filter_variants(input_vcf, output_vcf):
    """Filter variants based on quality and depth."""
    command = (
        f"bcftools filter -i 'QUAL>20 && DP>10' "
        f"{input_vcf} -Oz -o {output_vcf}"
    )
    run_command(command)

def main():
    reference = sys.argv[1]
    normal_bam = sys.argv[2]
    cancer_bam = sys.argv[3]
    # Normal tissue processing
    #normal_vcf = "normal_tissue_variant_pipeline.vcf.gz"
    #normal_filtered_vcf = "normal_tissue_variant_pipeline_filtered.vcf.gz"
    normal_vcf = normal_bam[:-4]+".vcf.gz"
    normal_filtered_vcf = normal_bam[:-4]+"_filtered.vcf.gz"
    call_variants(normal_bam, reference, normal_vcf)
    filter_variants(normal_vcf, normal_filtered_vcf)

    # Cancer tissue processing
    #cancer_vcf = "cancer_tissue_variant_pipeline.vcf.gz"
    #cancer_filtered_vcf = "cancer_tissue_variant_pipeline_filtered.vcf.gz"
    cancer_vcf = cancer_bam[:-4]+".vcf.gz"
    cancer_filtered_vcf = cancer_bam[:-4]+"_filtered.vcf.gz"
    call_variants(cancer_bam, reference, cancer_vcf)
    filter_variants(cancer_vcf, cancer_filtered_vcf)

if __name__ == "__main__":
    main()
