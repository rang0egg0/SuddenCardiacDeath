#!/bin/bash -l
#SBATCH --job-name=plinkQC.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=40g 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=stock526@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e plinkQC.%j.e 
#SBATCH -o plinkQC.%j.o 

# Load PLINK module version 1.90b6.10
module load plink/1.90b6.10

# Change directory to where your data is located
cd /home/durwa004/stock526/Aim3/case_control/data

# Perform QC without HWE (Hardy-Weinberg Equilibrium)
plink \
 --chr-set 31 --nonfounders --allow-no-sex --allow-extra-chr \
 --file SCD_CaseControl.pheno \  # Input file name
 --maf 0.01 --geno 0.1 \  # Minimum allele frequency and genotype missingness thresholds
 --make-bed --out SCD_CaseControl.pheno.QC  # Output file name prefix
