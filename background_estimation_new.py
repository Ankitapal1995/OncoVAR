import statistics
import sys

def parse_vcf(vcf_file):
    """Parse a VCF file to extract necessary details."""
    dp_values = []
    af_values = []
    records = []

    with open(vcf_file, "r") as f:
        for line in f:
            line = line.strip("\n")
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    header = line.split("\t")
                continue
            else:
                # Extract DP from INFO field
                info_field = line.split("\t")[7]
                dp = int(info_field.split(";")[1][3:])

                # Extract AF from FORMAT field
                format_field = line.split("\t")[9]
                if len(format_field.split(":")[2].split(",")) == 1:
                    af = float(format_field.split(":")[2])
                else:
                    af = float(format_field.split(":")[2].split(",")[1])
                # Append extracted values
                dp_values.append(dp)
                af_values.append(af)

                # Store the mutation record
                chrom = line.split("\t")[0]
                pos = line.split("\t")[1]
                ref = line.split("\t")[3]
                alt = line.split("\t")[4]
                records.append([chrom, pos, ".", ref, alt, dp, af])

    return dp_values, af_values, records

def write_summary(output_file, total_bases, dp_values, af_values):
    """Generate a summary report."""
    total_mutations = len(dp_values)
    total_reads=sum(dp_values)
    mutation_rate = total_mutations / total_bases if total_bases > 0 else 0
    median_background_mutation = statistics.median(af_values) if af_values else 0
    median_read_depth = statistics.median(dp_values) if dp_values else 0
    mean_alt_allele_fraction = statistics.mean(af_values) if af_values else 0
    stdev_alt_allele_fraction = statistics.stdev(af_values) if len(af_values) > 1 else 0
    reads_per_million = total_reads / 1000000 

    with open(output_file, "w") as f:
        f.write(f"Total mutations: {total_mutations}\n")
        f.write(f"Mutation rate: {mutation_rate}\n")
        f.write(f"Median background mutation level: {median_background_mutation}\n")
        f.write(f"Median read depth: {median_read_depth}\n")
        f.write(f"Mean alt allele fraction: {mean_alt_allele_fraction}\n")
        f.write(f"Standard deviation of alt allele fraction: {stdev_alt_allele_fraction}\n")
        f.write(f"Reads per million required to confidently call a mutation: {reads_per_million}\n")


def write_detailed_output(output_file, records):
    """Write detailed output for each mutation."""
    with open(output_file, "w") as f:
        f.write("#CHROM\tPOS\tID\tREF\tALT\tDP\tAF\n")
        for record in records:
            f.write("\t".join(map(str, record)) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python background_estimation.py <vcf_file> <total_bases> <output_prefix>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    total_bases = int(sys.argv[2])
    output_prefix = sys.argv[3]

    dp_values, af_values, records = parse_vcf(vcf_file)
    write_summary(f"{output_prefix}_summary.txt", total_bases, dp_values, af_values)
    write_detailed_output(f"{output_prefix}_detailed.txt", records)

