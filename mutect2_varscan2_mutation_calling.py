import os
import subprocess
import sys
import gzip
def run_command(command, description):
    """Run a shell command and print its description."""
    print(f"\nRunning: {description}\n{' '.join(command)}")
    subprocess.run(command, check=True)

def createRefdictionary(reference):
    "creating dictionary for reference"
    run_command(["/data/ankita/PB/project2/gatk-4.6.1.0/gatk", "CreateSequenceDictionary",
                "-R",reference],"Generating reference dictionary")
    
def generate_normal_sample_vcf(reference, normal_bam, normal_sample_vcf):
     run_command([
        "/data/ankita/PB/project2/gatk-4.6.1.0/gatk","Mutect2",
        "-R", reference,
        "-I", normal_bam,
        "-O", normal_sample_vcf
        ], "Generating normal sample vcf")

def run_mutect2(reference, tumor_bam, normal_bam, normal_sample, mutect_vcf, normal_sample_vcf,):
    """Run GATK Mutect2 for somatic variant calling."""
    #normal_sample = os.path.splitext(os.path.basename(normal_bam))[0]  # Remove '.bam' extension
    run_command([
        "/data/ankita/PB/project2/gatk-4.6.1.0/gatk", "Mutect2",
        "-R", reference,
        "-I", tumor_bam,
        "-I", normal_bam,
        "--normal", normal_sample,
        "--panel-of-normals", normal_sample_vcf,
        "-O", mutect_vcf
    ], "Running Mutect2")

def filter_mutect2(reference,mutect_vcf, filtered_mutect_vcf):
    """Filter low-confidence variants from Mutect2 output."""
    run_command([
        "/data/ankita/PB/project2/gatk-4.6.1.0/gatk", "FilterMutectCalls",
        "-R", reference,
        "-V", mutect_vcf,
        "-O", filtered_mutect_vcf
    ], "Filtering Mutect2 calls")


def generate_combined_pileup(reference, normal_bam, tumor_bam, combined_pileup_file):
    run_command([
        "samtools", "mpileup",
        "-f", reference,
        "-q", "20",     
        "-B",
        "-o", combined_pileup_file,
        normal_bam, tumor_bam
    ], "Generating combined mpileup with Samtools")


def run_varscan2(combined_pileup_file, varscan_vcf_prefix):
    """Run VarScan2 to identify somatic mutations from a combined mpileup."""
    run_command([
        "varscan", "somatic",
        combined_pileup_file, varscan_vcf_prefix,
        "--mpileup", "1",
        "--min-coverage", "30",
        "--min-var-freq", "0.05",
        "--somatic-p-value", "0.05",
        "--output-vcf 1"
    ], "Running VarScan2 somatic mutation calling")

def compressed_vcf_outputs(varscan_vcf_prefix_snp,varscan_vcf_prefix_indel):
    "compressed varsan snp and indels"
    run_command(["bgzip",varscan_vcf_prefix_snp], f"Compressing varscan snp vcf,")
    run_command(["bgzip",varscan_vcf_prefix_indel], f"Compressing varscan indel vcf,")
def indexing_vcfs(comp_varscan_vcf_prefix_snp,comp_varscan_vcf_prefix_indel):
    run_command(["bcftools", "index", comp_varscan_vcf_prefix_snp], f"indexing varscan snp vcf")
    run_command(["bcftools", "index", comp_varscan_vcf_prefix_indel], f"indexing varscan indel vcf")



def merge_vcfs(varscan_vcf_prefix_snp, varscan_vcf_prefix_indel, varscan_merged_vcf):
    with open(varscan_merged_vcf, 'w') as outfile:
        # Read all lines from the first VCF file
        with open(varscan_vcf_prefix_snp, 'rt') as snp_file:
            outfile.writelines(snp_file.readlines())
        # Read non-header lines from the second VCF file
        with open(varscan_vcf_prefix_indel, 'rt') as indel_file:
            for line in indel_file:
                if not line.startswith("#"):
                    outfile.write(line)


def estimate_background_mutation(combined_vcf, background_file):
    """Perform background mutation estimation."""
    # Placeholder for a user-defined estimation script or logic
    run_command([
        "python", "background_estimation.py", combined_vcf, background_file
    ], "Estimating background mutations")

def main():
    if len(sys.argv) != 5:
        print("Usage: python code.py reference.fasta tumor.bam normal.bam normal_sample")
        sys.exit(1)

    reference = sys.argv[1]
    tumor_bam = sys.argv[2]
    normal_bam = sys.argv[3]
    normal_sample = sys.argv[4]

    # Define intermediate and output file names
    normal_sample_vcf="normal_sample.vcf"
    pon_file="pon.vcf.gz"
    mutect_vcf = "mutect2_output.vcf"
    filtered_mutect_vcf = "filtered_mutect2_output.vcf"
    combined_pileup_file ="tumor_normal.mpileup"
    varscan_vcf_prefix = "varscan_output"
    varscan_vcf = f"{varscan_vcf_prefix}.somatic.vcf"
    varscan_vcf_prefix_snp = f"{varscan_vcf_prefix}.snp.vcf"
    varscan_vcf_prefix_indel = f"{varscan_vcf_prefix}.indel.vcf"
    comp_varscan_vcf_prefix_snp = varscan_vcf_prefix+".snp.vcf.gz"
    comp_varscan_vcf_prefix_indel = varscan_vcf_prefix+".indel.vcf.gz"
    comp_varscan_vcf_prefix_snp_indexed = varscan_vcf_prefix+".snp.vcf.gz.csi"
    comp_varscan_vcf_prefix_indel_indexed = varscan_vcf_prefix+".indel.vcf.gz.csi"
    varscan_merged_vcf ="varscan_final_output.vcf"
    compressed_mutect2="filtered_mutect2_output.vcf.gz"
    compressed_varscan= "varscan_final_output.vcf.gz"
    compressed_mutect2_indexed = compressed_mutect2+".csi"
    compressed_varscan_indexed = compressed_varscan+".csi"
    combined_vcf = "combined_somatic_mutations.vcf"
    background_file = "background_mutation_estimation.txt"
    
   # Running mutect2 pipeline

    createRefdictionary(reference)

    generate_normal_sample_vcf(reference, normal_bam, normal_sample_vcf)

    run_mutect2(reference, tumor_bam, normal_bam, normal_sample, mutect_vcf, normal_sample_vcf)

    filter_mutect2(reference,mutect_vcf, filtered_mutect_vcf)

   # Running varscan2 pipeline
   
    generate_combined_pileup(reference, normal_bam, tumor_bam, combined_pileup_file)

    run_varscan2(combined_pileup_file, varscan_vcf_prefix)

    merge_vcfs(varscan_vcf_prefix_snp, varscan_vcf_prefix_indel, varscan_merged_vcf)

    merge_vcfs(filtered_mutect_vcf, varscan_merged_vcf, combined_vcf)
    
    #estimate_background_mutation(combined_vcf, background_file)
    print(f"Pipeline completed. Final combined VCF: {combined_vcf}")
    print(f"Background mutation estimation file: {background_file}")

if __name__ == "__main__":
    main()

