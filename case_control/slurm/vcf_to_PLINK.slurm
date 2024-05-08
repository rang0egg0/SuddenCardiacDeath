#!/bin/bash -l
#SBATCH --job-name=plink.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=40g 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=stock526@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e plink.%j.e 
#SBATCH -o plink.%j.o 

module load plink/1.90b6.10 
cd /home/durwa004/stock526/Aim3/case_control/data

#Convert VCF to plink binary files
plink\
 --chr-set 31 --nonfounders --double-id --allow-no-sex\
 --vcf SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf.gz\
 --make-bed --allow-extra-chr --vcf-half-call missing\
 --out SCD_CaseControl 

#Add in phenotype data
plink\
 --chr-set 31 --nonfounders --allow-no-sex --allow-extra-chr\
 --bfile SCD_CaseControl\
 --pheno SCD_CC_TFAM.tfam --pheno-name PHENOTYPE --recode\
 --out SCD_CaseControl.pheno
