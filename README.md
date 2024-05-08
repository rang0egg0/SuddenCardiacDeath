# SuddenCardiacDeath

# CANDIDATE GENE DISCOVERY 

### Step 0.1 Mine candidate genes from keyword-based databases
Keywords used: "sudden cardiac death", "cardiac arrhythmia", and "sudden arrhythmic death" 
Programs used: Phenolyzer (https://phenolyzer.wglab.org/), ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/), OpenTargets (https://www.opentargets.org/), OMIM (https://www.omim.org/)
Download outputs on local computer and then upload them to MSI. 

### Step 0.2 (OpenOnDemand, R)Candidate gene data wrangling and finding Union 
Candidate_Genes/jan24can_gen.R
#### Outputs two files Narrow (present in all programs) and Broad (present in two or more programs) candidate genes 



# CASE CONTROL WORK  

### Step 1.1: (slurm, bcftools) Subset cases and controls from population VCF
sbatch case_control/slurm/subsetting_population_by_CasesControls.slurm


### Step 1.2: (slurm, SnpSift) CaseControl annotations
sbatch case_control/slurm/SnpSift_case_control.slurm
#### calculates four p-values for each variant based on four different models of genetic expression (dominant, recessive, codominant/genotypic, and Cochran Armitage). P-values are appended to the end of the INFO field


### Step 1.3: (slurm, PLINK) Create binary PLINK files for data
sbatch case_control/slurm/vcf_to_PLINK.slurm
#### creates PLINK binary outputs from VCF and adds in phenotype (case/control status) data


### Step 1.4: (slurm, PLINK) Perform quality control on PLINK output
sbatch case_control/slurm/plinkQC.slurm
#### performs QC without regard to Hardy-Weinberg equilibrium and removes variants that do not meet QC standards


### Step 1.5: (slurm, PLINK) Prune variants in linkage disequilibrium 
sbatch case_control/slurm/plink_prune.slurm
#### calculates which variants are in LD and then creates a new file that includes only variants that are not in LD


### Step 1.6: (Python) Get type of variant
sbatch case_control/python/get_type_of_variant.py


### Step 1.7: (Python) Prune variants that do not meet P-Value thresholds 
sbatch case_control/python/SnpSiftResults.py
#### Retains only variants where at least one of the P-Values calculated in Step 1.2 meets the significance threshold. The significance threshold is based on a Bonferroni correction (0.05/x), x = the total number of variants. The value for "x" can be found in the .log file for PLINK analysis in Step __. 


### Step 1.8 (

## GWAS WORK 

### Step 2.1: (slurm, PLINK) Determine associations based on PLINK analysis
sbatch case_control/GWAS/plink_association.slurm
#### Uses binary output files from Step 1.4 as input


### Step 2.2 (slurm, GEMMA) Create a relationship matrix to be used in GWAS
sbatch case_control/GWAS/gemma_gwas.slurm
#### Uses binary output files from Step 1.5 as input 


### Step 2.3 (OpenOnDemand, R) Create the Manhattan Plots
case_control/GWAS/manhattan_plots.R
#### Update file paths and significance thresholds in R file. Uses binary outputs from Step 2.1 and Step 2.2 as inputs. The significance threshold for et_sig is based on the number of variants after LD correction, found in the .log file for Step 1.__. The significance threshold for gw_sig is based on the total number of variants, which can be found in the .log file for Step 1.__. 
#### It is best to run this R script (and any R script) on the MSI OpenOnDemand platform (https://ondemand.msi.umn.edu/pun/sys/dashboard/) 







# VARIANT ANALYSIS BY CANDIDATE GENES

