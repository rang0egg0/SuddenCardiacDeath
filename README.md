# SuddenCardiacDeath

## CASE CONTROL WORK # 

### Step 1: (bcftools) Subset cases and controls from population VCF
sbatch case_control/slurm/subsetting_population_by_CasesControls.slurm


### Step 2: (SnpSift) CaseControl annotations
sbatch case_control/slurm/SnpSift_case_control.slurm
#### calculates four p-values for each variant based on four different models of genetic expression (dominant, recessive, codominant/genotypic, and Cochran Armitage). P-values are appended to the end of the INFO field


### Step 3: (PLINK) Create binary PLINK files for data
sbatch case_control/slurm/vcf_to_PLINK.slurm
#### creates PLINK binary outputs from VCF and adds in phenotype (case/control status) data


### Step 4: (PLINK) Perform quality control on PLINK output
sbatch case_control/slurm/plinkQC.slurm
#### performs QC without regards to Hardy-Weinberg equilibrium and removes variants that do not meet QC standards


### Step 5: (PLINK) Prune variants in linkage disequalibrium 
sbatch case_control/slurm/plink_prune.slurm
