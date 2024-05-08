#!/bin/bash -l
#SBATCH --job-name=SnpEffCC.slurm
#SBATCH --nodes=1 # Specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 # Not sure if this should be 2 or 1 for consistency
#SBATCH --mem=10g 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=stock526@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e SnpEffCC.%j.e 
#SBATCH -o SnpEffCC.%j.o 

# Activate the Conda environment named joysmamba
conda activate joysmamba

# Unzip the input VCF file
gunzip /home/durwa004/stock526/Aim3/case_control/data/SCD_CaseControl.goldenPath.vcf.gz

# Load the Java OpenJDK module (required to run SnpSift)
module load java/openjdk-17.0.2

# Move to the snpEff directory
cd /home/durwa004/durwa004/packages/snpEff/

# Run SnpSift.jar to perform case-control annotation
java -jar SnpSift.jar caseControl \
    /home/durwa004/stock526/Aim3/case_control/data/SCD_CaseControl.goldenPath.vcf \
    -tfam /home/durwa004/stock526/Aim3/case_control/data/SCD_CC_TFAM.tfam \
    > /home/durwa004/stock526/Aim3/case_control/data/SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf

# Move back to the data directory
cd /home/durwa004/stock526/Aim3/case_control/data

# Compress the annotated VCF file using bgzip
bgzip SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf
bgzip SCD_CaseControl.goldenPath.vcf

# Create a tabix index file for the annotated VCF
tabix SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf.gz

