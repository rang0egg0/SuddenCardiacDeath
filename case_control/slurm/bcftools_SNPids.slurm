#!/bin/bash
#SBATCH --job-name=vcf_id_assignment
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stock526@umn.edu

conda activate joysmamba

# Load the necessary modules
module load bcftools

# Define input and output file paths
INPUT_VCF="SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf.gz"
OUTPUT_VCF="SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.ID.vcf.gz"

# Add a unique ID to each variant in the VCF file based on chromosome and position
bcftools annotate --set-id +'%CHROM\_%POS\_%ALT' $INPUT_VCF -Oz -o $OUTPUT_VCF
