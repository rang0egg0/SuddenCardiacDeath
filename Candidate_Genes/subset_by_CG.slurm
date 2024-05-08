#!/bin/bash -l
#SBATCH --job-name=cg_comp.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=60g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stock526@umn.edu
#SBATCH -A durwa004
#SBATCH -e cg_comp2.%j.e
#SBATCH -o cg_comp2.%j.o
#SBATCH --cpus-per-task=1
#SBATCH -p msismall,msilarge



# Activate Conda environment and load bcftools module
conda activate joysmamba
module load bcftools

# Change to the directory containing the VCF files
cd /home/durwa004/stock526/Aim3/candidate_gene/data/

# Filter cases VCF by broad candidate genes
bcftools view cases.goldenPath.vcf.gz \
 -R /home/durwa004/stock526/CandidateGenes/results/broad_candidate_genes.txt \
 -Oz -o cases_broadCG.goldenPath.vcf.gz

# Filter cases VCF by narrow candidate genes
bcftools view cases.goldenPath.vcf.gz \
 -R /home/durwa004/stock526/CandidateGenes/results/narrow_candidate_genes.txt \
 -Oz -o cases_narrowCG.goldenPath.vcf.gz

# Filter genpop VCF by broad candidate genes
bcftools view genpop.goldenPath.vcf.gz \
 -R /home/durwa004/stock526/CandidateGenes/results/broad_candidate_genes.txt \
 -Oz -o genpop_broadCG.goldenPath.vcf.gz

# Filter genpop VCF by narrow candidate genes
bcftools view genpop.goldenPath.vcf.gz \
 -R /home/durwa004/stock526/CandidateGenes/results/narrow_candidate_genes.txt \
 -Oz -o genpop_narrowCG.goldenPath.vcf.gz

# Index filtered VCF files using tabix
tabix genpop_narrowCG.goldenPath.vcf.gz
tabix genpop_broadCG.goldenPath.vcf.gz
tabix cases_narrowCG.goldenPath.vcf.gz
tabix cases_broadCG.goldenPath.vcf.gz
