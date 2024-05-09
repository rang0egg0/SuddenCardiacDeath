#!/bin/bash -l
#SBATCH --job-name=vcftoplink.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=40g 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=stock526@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e vcftoplink.%j.e 
#SBATCH -o vcftoplink.%j.o 

# Load the PLINK software module version 1.90b6.10
module load plink/1.90b6.10

# Change directory to where your data is located
cd /home/durwa004/stock526/Aim3/case_control/data

# Run PLINK with the following options:
# --chr-set 31: Set the chromosome code to 31
# --nonfounders: Include non-founders in the analysis
# --double-id: Allow duplicate IDs
# --allow-no-sex: Allow samples with missing sex information
# --vcf: Input VCF file name
# --make-bed: Convert VCF to PLINK binary format
# --allow-extra-chr: Allow extra chromosomes in the data
# --vcf-half-call missing: Treat half-calls in VCF files as missing data
# --pheno: Input phenotype file name
# --pheno-name PHENOTYPE: Specify the column name for phenotype data
# --recode: Recode the data after adding phenotype information
# --out: Output file name prefix

plink \
 --chr-set 31 --nonfounders --double-id --allow-no-sex \
 --vcf SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf.gz \
 --make-bed --allow-extra-chr --vcf-half-call missing \
 --pheno SCD_CC_TFAM.tfam --pheno-name PHENOTYPE --recode \
 --out SCD

