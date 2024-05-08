import argparse  # Importing the argparse module for command-line argument parsing
import os  # Importing the os module for file and path operations
import gzip  # Importing the gzip module for reading gzipped files

def make_arg_parser():  # Function to create an argument parser for the script
    parser = argparse.ArgumentParser(
            prog="GeneratePBS.py",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            "-d", "--data",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="input vcf.gz file [required]")  # Argument for input VCF.gz file
    parser.add_argument(
            "-p", "--prog",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Variant annotator: SnpEff or VEP [required]")  # Argument for variant annotator
    return parser

if __name__ == '__main__':
    parser = make_arg_parser()  # Create argument parser object
    args = parser.parse_args()  # Parse command-line arguments

    data = os.path.abspath(args.data)  # Get absolute path of input VCF.gz file
    prog = args.prog  # Get selected variant annotator (SnpEff or VEP)
    if prog == "SnpEff":
        prog2 = "ANN="  # Annotation prefix for SnpEff
    elif prog == "VEP":
        prog2 = "CSQ="  # Annotation prefix for VEP

# Goal is to convert the snpeff output to a usable text file for analysis, combining with annovar output.
with gzip.open(data, "rt") as input_file, open(f"{prog}.hml.txt", "w") as output_file, open(f"{prog}.genes.txt", "w") as output2:
    high = 0  # Counter for high impact variants
    mod = 0  # Counter for moderate impact variants
    low = 0  # Counter for low impact variants
    for line in input_file:  # Iterate through each line in the input VCF.gz file
        line = line.rstrip("\n").split("\t")  # Split the line by tabs

        # Process header lines
        if "#" in line[0]:
            if prog == "VEP":
                if "##INFO=<ID=CSQ" in line[0]:
                    line = line[0].split("Allele|")
                    line = line[1].split("|")
                    # Write header line to VEP output file
                    print("chrom\tpos\tref\talt\tAC\tAF",  "\t".join(line[0:30]), sep="\t", file=output_file)
                else:
                    continue  # Skip other header lines
            elif prog == "SnpEff":
                if "##INFO=<ID=ANN" in line[0]:
                    line = line[0].split("Allele | ")
                    line = line[1].split(" | ")
                    for i, v in enumerate(line):
                        if "/" in v:
                            v = v.split(" / ")
                            line[i] = v[0] + "_" + v[1]
                    # Write header line to SnpEff output file
                    print("chrom\tpos\tref\talt\tAC\tAF", "\t".join(line), "LOF\tgene_name\tgene_ID\tno_transcripts\tpercent_transcripts_affected\tcases_het\tcases_hom\tcases_count\tcontrols_het\tcontrols_hom\tcontrols_count\tCC_trend\tCC_geno\tCC_all\tCC_dom\tCC_rec", sep="\t", file=output_file)
                else:
                    continue  # Skip other header lines
        else:  # Process variant lines
            chr_pos = line[0] + ":" + line[1] + ":" + line[3]  # Create chr:pos:ref string
            al = line[4].split(",")  # Split alternative alleles by comma
            ab = line[7].split(";")  # Split INFO field by semicolon
            bc = ab[0].split("AC=")[1]  # Extract allele count (AC)
            cd = ab[1].split("AF=")[1]  # Extract allele frequency (AF)
            ef = line[7].split(prog2)  # Split INFO field by annotation prefix
            if len(ef) > 1:
                ef = ef[1].split(";")
                ef = ef[0].split(",")
                fg = ef[0].split("|")
                plant = []
                for item, value in enumerate(fg[1:]):
                    if item == 30:
                        break
                    elif item == 3 or item == 2 or item == 5:
                        value = value.split(".")
                        plant.append(value[0])
                        print(value[0], file=output2)
                    elif value == "":
                        plant.append("N")
                    else:
                        plant.append(value)
                # Check impact and update impact counters
                if fg[2] == "HIGH":
                    high += 1
                elif fg[2] == "MODERATE":
                    mod += 1
                elif fg[2] == "LOW":
                    low += 1
                chr_pos1 = chr_pos.split(":")
                aa = line[7].split("LOF=(")
                if prog == "SnpEff":
                    # Extract additional SnpEff-specific annotations
                    cases = line[7].split("Cases=")[1].split(";")[0].split(",")
                    controls = line[7].split("Controls=")[1].split(";")[0].split(",")
                    trend = line[7].split("CC_TREND=")[1]
                    geno = line[7].split("CC_GENO=")[1]
                    cc_all = line[7].split("CC_ALL=")[1]
                    dom = line[7].split("CC_DOM=")[1]
                    rec = line[7].split("CC_REC=")[1]
                    # Write variant information to SnpEff output file
                    if len(aa) > 1:
                        bb = aa[1].split(")")[0].split("|")
                        LOF = "yes"
                        print("\t".join(chr_pos1), line[4], bc, cd, "\t".join(plant), LOF, "\t".join(bb), "\t".join(cases), "\t".join(controls), dom, rec, cc_all, geno, trend, sep="\t", file=output_file)
                    else:
                        LOF = "no"
                        NL = list("N" * 4)
                        print("\t".join(chr_pos1), line[4], bc, cd, "\t".join(plant), LOF, "\t".join(NL), "\t".join(cases), "\t".join(controls), dom, rec, cc_all, geno, trend, sep="\t", file=output_file)
                elif prog == "VEP":
                    # Write variant information to VEP output file
                    print("\t".join(chr_pos1), line[4], bc, cd, "\t".join(plant), sep="\t", file=output_file)

    # Print summary of impact variant counts
    print(f"High impact variants: {high}\nModerate impact variants: {mod}\nLow impact variants: {low}")
