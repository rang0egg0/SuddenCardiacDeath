#!/bin/bash -l
#SBATCH --job-name=gemma_geas.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=10g 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=stock526@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e gwas_e.%j.e 
#SBATCH -o gwas_o.%j.o 

cd /home/durwa004/stock526/Aim3/GWAS/  # Change directory to the GWAS directory

# Calculate the genetic relationship matrix (GRM) using GEMMA
/home/durwa004/durwa004/packages/gemma \
 -bfile /home/durwa004/stock526/Aim3/case_control/data/SCD_CaseControl.pheno.QC.LD.pruned.final \  # Input PLINK binary file for GRM calculation
 -gk 2 \  # Use the VanRaden method for GRM calculation
 -o SCD_CaseControl.pheno.QC.LD.pruned.final.GRM  # Output prefix for GRM files

# Run genome-wide association analysis using GEMMA
/home/durwa004/durwa004/packages/gemma \
 -bfile /home/durwa004/stock526/Aim3/case_control/data/SCD_CaseControl.pheno.QC.LD.pruned.final \  # Input PLINK binary file for association analysis
 -k ./output/SCD_CaseControl.pheno.QC.LD.pruned.final.GRM.sXX.txt \  # Input GRM file (sXX.txt format)
 -fa 4 \  # Specify the number of threads for analysis
 -o gemma_gwas_SCD_CaseControl  # Output prefix for association analysis results
