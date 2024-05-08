#!/bin/bash -l
#SBATCH --job-name=cg_comp.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=60g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stock526@umn.edu
#SBATCH -A durwa004
#SBATCH -e cg_comp.%j.e
#SBATCH -o cg_comp.%j.o
#SBATCH --cpus-per-task=1
#SBATCH -p msismall,msilarge

conda activate joysmamba  # Activate Conda environment named "joysmamba"
module load bcftools  # Load the bcftools module

# Subset by cases
bcftools view ~/joint_call_passing.annotated.goldenPath.20230726.vcf.gz \
 -S /home/durwa004/stock526/Aim3/candidate_gene/data/SCD_case_only_IDs.txt \
 -Oz -o /home/durwa004/stock526/Aim3/candidate_gene/data/cases.goldenPath.vcf.gz

# Subset gen_pop so that it does not include cases
bcftools view ~/joint_call_passing.annotated.goldenPath.20230726.vcf.gz \
 -S ^/home/durwa004/stock526/Aim3/candidate_gene/data/SCD_case_only_IDs.txt \
 -Oz -o /home/durwa004/stock526/Aim3/candidate_gene/data/genpop.goldenPath.vcf.gz

# Compress output files
cd /home/durwa004/stock526/Aim3/candidate_gene/data/  # Change to the data directory
tabix genpop.goldenPath.vcf.gz  # Index the genpop VCF file
tabix cases.goldenPath.vcf.gz  # Index the cases VCF file
