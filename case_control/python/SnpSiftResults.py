vcf = "ANSAD.SnpEff.hml.txt"  # Input VCF file name

# Significance threshold based on Bonferroni correction: 0.05 divided by the total number of variants
# The denominator should be updated based on the specific dataset
significance = 0.05 / 8682578  # Example significance threshold

count = 0  # Counter for variants passing the significance threshold

# Open the input VCF file for reading and the output file for writing
with open(vcf, "r") as input_f, open("ANSAD_SnpSift_case_control.txt", "w") as output_f:
    for line in input_f:  # Iterate through each line in the input VCF file
        line = line.rstrip("\n").split("\t")  # Split the line by tabs

        if "chrom" == line[0]:  # Check if the line is the header line
            print("\t".join(line), file=output_f)  # Write the header line to the output file
        else:  # Process variant lines
            # Check if any of the significance values in the last five columns are below the significance threshold
            if float(line[-1]) < significance or float(line[-2]) < significance or float(line[-3]) < significance or float(line[-4]) < significance or float(line[-5]) < significance:
                print("\t".join(line), file=output_f)  # Write the variant line to the output file
                count += 1  # Increment the count of variants passing the threshold

# Print the total count of variants passing the significance threshold
print(f"Total variants passing significance threshold: {count}")
