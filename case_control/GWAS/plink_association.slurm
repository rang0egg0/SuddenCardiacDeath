#!/bin/bash -l
#SBATCH --job-name=plink_assocation.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=10g 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=stock526@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e plink_ass.%j.e 
#SBATCH -o plink_ass.%j.o 

# Load PLINK module version 1.90b6.10
module load plink/1.90b6.10

# Change directory to where your data is located
cd /home/durwa004/stock526/Aim3/case_control/data

# Perform association analysis using PLINK
plink \
 --chr-set 31 --nonfounders --allow-no-sex --allow-extra-chr \
 --bfile SCD_CaseControl.pheno.QC \  # Input PLINK binary file name
 --assoc \  # Perform association analysis
 --out SCD_CaseControl.pheno.QC.assoc  # Output file name prefix for association results
