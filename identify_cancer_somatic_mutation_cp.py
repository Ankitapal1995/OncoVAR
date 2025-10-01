import vcfpy
import csv
from sys import argv
_,normal_vcf,cancer_vcf=argv
def parse_vcf(vcf_file):
    """
    Parse the VCF file and return a dictionary of variants with additional fields 
    (chromosome, position, reference, alternate, DP, AC, AN).
    """
    variants = {}
    reader = vcfpy.Reader.from_path(vcf_file)
    for record in reader:
        # Extract required fields
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = tuple([alt.value for alt in record.ALT])
        dp = record.INFO.get("DP", "NA")  # Get DP (read depth), default to "NA" if not present
        ac = record.INFO.get("AC", [0])  # Get AC (allele count), default to [0] if not present
        an = record.INFO.get("AN", 1)    # Get AN (allele number), default to 1 if not present
        
        # Store variant details in a dictionary
        variants[(chrom, pos, ref, alt)] = {
            "DP": dp,
            "AC": ac[0],  # Assuming AC is a single value
            "AN": an
        }
    return variants

def compare_variants(normal_vcf, cancer_vcf):
    """
    Compare variants between normal and cancer samples and find somatic variants.
    Calculate AF (Allele Frequency) for somatic variants.
    """
    normal_variants = parse_vcf(normal_vcf)
    cancer_variants = parse_vcf(cancer_vcf)

    somatic_variants = {}
    for variant, details in cancer_variants.items():
        if variant not in normal_variants:
            # Calculate AF as AC / AN
            ac = details["AC"]
            an = details["AN"]
            af = round(ac / an, 4) if an > 0 else "NA"
            somatic_variants[variant] = {
                "DP": details["DP"],
                "AC": ac,
                "AN": an,
                "AF": af
            }
    return somatic_variants


# Compare and find somatic variants
somatic_variants = compare_variants(normal_vcf, cancer_vcf)

# Write results to a CSV file
output_csv = "output.csv"
with open(output_csv, mode="w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    # Write header
    writer.writerow(["Chromosome", "Position", "Ref", "Alt", "DP", "AC", "AN", "AF"])
    # Write data
    for variant, details in somatic_variants.items():
        chrom, pos, ref, alt = variant
        writer.writerow([chrom, pos, ref, ", ".join(alt), details["DP"], details["AC"], details["AN"], details["AF"]])

# Print results to the console
for variant, details in somatic_variants.items():
    chrom, pos, ref, alt = variant
    print(f"Somatic Variant: Chromosome: {chrom}, Position: {pos}, Ref: {ref}, Alt: {alt}, "
          f"DP: {details['DP']}, AC: {details['AC']}, AN: {details['AN']}, AF: {details['AF']}")

print(f"\nSomatic variant results saved to {output_csv}.")

