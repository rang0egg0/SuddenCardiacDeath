#!/bin/bash -l
#SBATCH --job-name=subset_cases.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=60g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stock526@umn.edu
#SBATCH -A durwa004
#SBATCH -e 2sub_cases.%j.e
#SBATCH -o 2sub_cases%j.o
#SBATCH --cpus-per-task=1
#SBATCH -p msismall,msilarge

# Activate conda environment (can use your own or someone else's)
conda activate joysmamba

# Load BCFtools (does not require any sort of external downloading unlike other packages like GEMMA for example)
module load bcftools

# Subsetting from decomposed population VCF that contains VEP and SnpEff annotations and PhyloP scores
bcftools view /scratch.global/marlo072/Joy/outputs/final/joint_call_passing.annotated.goldenPath.20230726.vcf.gz \
# Subset by case and control IDs
 -S /home/durwa004/stock526/Aim3/case_control/data/cases_AND_controls_IDs.txt \
# Output as a bgzipped VCF
 -Oz -o /home/durwa004/stock526/Aim3/case_control/data/2SCD_CaseControl2.goldenPath.vcf.gz

# Move into directory with output
cd  /home/durwa004/stock526/Aim3/case_control/data

# Create a tabix index file for output
tabix SCD_CaseControl.goldenPath2.vcf.gz
