#!/bin/bash -l
#SBATCH --job-name=plink_LD.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=10g 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=stock526@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e plink_LD.%j.e 
#SBATCH -o plink_LD.%j.o 


# Load PLINK module version 1.90b6.10
module load plink/1.90b6.10

# Change directory to where your data is located
cd /home/durwa004/stock526/Aim3/case_control/data/

# Perform LD pruning using PLINK
plink \
 --chr-set 31 --nonfounders --allow-no-sex --allow-extra-chr \
 --bfile SCD_CaseControl.pheno.QC \  # Input file name in PLINK binary format
 --indep-pairwise 100 25 0.2 \  # LD pruning parameters (window size, step, and r^2 threshold)
 --out SCD_CaseControl.pheno.QC.LD.pruned  # Output file name prefix for LD pruned data

# Keep only SNPs not in LD based on LD-pruned data
plink \
 --chr-set 31 --nonfounders --allow-no-sex --allow-extra-chr \
 --bfile SCD_CaseControl.pheno.QC \  # Input file name in PLINK binary format
 --extract SCD_CaseControl.pheno.QC.LD.pruned.prune.in \  # List of SNPs not in LD from the LD-pruned data
 --make-bed --out SCD_CaseControl.pheno.QC.LD.pruned.final  # Output file name prefix for final pruned data
