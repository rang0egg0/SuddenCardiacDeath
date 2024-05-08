#!/bin/bash -l
#SBATCH --job-name=SnpEffCC.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=10g 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=stock526@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e SnpEffCC.%j.e 
#SBATCH -o SnpEffCC.%j.o 

conda activate joysmamba

gunzip /home/durwa004/stock526/Aim3/case_control/data/SCD_CaseControl.goldenPath.vcf.gz

module load java/openjdk-17.0.2
cd /home/durwa004/durwa004/packages/snpEff/

java -jar SnpSift.jar caseControl /home/durwa004/stock526/Aim3/case_control/data/SCD_CaseControl.goldenPath.vcf -tfam /home/durwa004/stock526/Aim3/case_control/data/SCD_CC_TFAM.tfam  > /home/durwa004/stock526/Aim3/case_control/data/SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf


cd /home/durwa004/stock526/Aim3/case_control/data

bgzip SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf
bgzip SCD_CaseControl.goldenPath.vcf

tabix SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf.gz
