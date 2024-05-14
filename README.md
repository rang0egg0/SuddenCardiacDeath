# SuddenCardiacDeath

# CANDIDATE GENE DISCOVERY 

### Step 0.1 Mine candidate genes from keyword-based databases
Keywords used: "sudden cardiac death", "cardiac arrhythmia", and "sudden arrhythmic death" 
Programs used: Phenolyzer (https://phenolyzer.wglab.org/), ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/), OpenTargets (https://www.opentargets.org/), OMIM (https://www.omim.org/)
Download outputs on local computer and then upload them to MSI. 

### Step 0.2 (OpenOnDemand, R) Candidate gene data wrangling and finding Union 
Outputs two files Narrow (present in all programs) and Broad (present in two or more programs) candidate genes 
```
jan24can_gen.R
```


# CASE CONTROL WORK  

### Step 1.1: (slurm, bcftools) Subset cases and controls from population VCF
#### Infile
joint_call_passing.annotated.goldenPath.20230726.vcf.gz
#### Outfile
SCD_CaseControl.goldenPath.vcf.gz 
```
sbatch subsetting_population_by_CasesControls.slurm
```

### Step 1.2: (slurm, SnpSift) CaseControl annotations
calculates four p-values for each variant based on four different models of genetic expression (dominant, recessive, codominant/genotypic, and Cochran Armitage). P-values are appended to the end of the INFO field
#### Infile
SCD_CaseControl.goldenPath.vcf.gz 
#### Outfile
SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf.gz
```
sbatch SnpSift_case_control.slurm
```

### Step 1.3 (slurm, BCFtools) Create unique IDs for SNPs
Annotates the "ID" column of the vcf to give each SNP a unique ID (chrom_pos_alt)
#### Infile
SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.vcf.gz
#### Outfile
SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.ID.vcf.gz
```
sbatch bcftools_generate_SNPids.slurm
tabix SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.ID.vcf.gz
```

### Step 1.4: (slurm, PLINK) Create binary PLINK files for data
creates PLINK binary outputs from VCF and adds in phenotype (case/control status) data
#### Infile
SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.ID.vcf.gz
#### Outfile
SCD.(bim, bed, fam, log, map, nosex, ped)
```
sbatch vcf_to_PLINK.slurm
```

### Step 1.5: (slurm, PLINK) Perform quality control on PLINK output
performs QC without regard to Hardy-Weinberg equilibrium and removes variants that do not meet QC standards
#### Infile
SCD.     
#### Outfile
SCD.QC.
```
cd PLINK
sbatch plinkQC.slurm
```

### Step 1.6.1: (slurm, PLINK) Prune variants in linkage disequilibrium 
calculates which variants are in LD putting the SNP_IDs in a .in and .out file
#### Infile
SCD.QC  
#### Outfile
SCD.QC.LDpruned.prune.in
SCD.QC.LDpruned.prune.out
SCD.QC.LDpruned.
```
sbatch plink_prune1.slurm
```
### Step 1.6.2: (slurm, PLINK) Prune variants in linkage disequilibrium 
Filters out variants that are in LD from SCD.QC. 
#### Infile
SCD.QC
SCD.QC.LDpruned.prune.in
#### Outfile
SCD.QC.LDpruned.final
```
sbatch plink_prune2.slurm
```

### Step 1.7 Plink back to VCF
#### Infile
SCD.QC.LDpruned.final

#### Outfile
SCD_final_genotypes.vcf.gz
```
sbatch plink_to_VCF.slurm
bgzip SCD_final_genotypes.vcf
tabix SCD_final_genotypes.vcf.gz
```

### Step 1.8 Obtain SNP IDs
#### Infile
SCD_final_genotypes.vcf.gz
#### Outfile
final_snp_ids.txt

### Step 1.8 Subset VCF by snpID
#### Infile
SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.ID.vcf.gz
#### Outfile
final.vcf.gz
```
bcftools view -i 'ID=@final_snp_ids.txt' SCD_CaseControl.goldenPath.decomposed.snpeff.snpsift.CC.ID.vcf.gz -o final.vcf
bgzip final.vcf
tabix final.vcf.gz

### Step 1.7: (Python) Get type of variant
```
ipython /home/durwa004/stock526/Aim3/case_control/scripts/python/get_type_of_variant.py -d /home/durwa004/stock526/Aim3/case_control/data/inal.vcf.gz -p SnpEff
```

### Step 1.8: (Python) Prune variants that do not meet P-Value thresholds 
Retains only variants where at least one of the P-Values calculated in Step 1.2 meets the significance threshold. The significance threshold is based on a Bonferroni correction (0.05/x), x = the total number of variants. The value for "x" can be found in the .log file for PLINK analysis in Step __. 
```
sbatch SnpSiftResults.py
```


## GWAS WORK 

### Step 2.1: (slurm, PLINK) Determine associations based on PLINK analysis
Uses binary output files from Step 1.4 as input
```
sbatch plink_association.slurm
```

### Step 2.2 (slurm, GEMMA) Create a relationship matrix to be used in GWAS
Uses binary output files from Step 1.5 as input 
```
sbatch gemma_gwas.slurm
```

### Step 2.3 (OpenOnDemand, R) Create the Manhattan Plots
Update file paths and significance thresholds in R file. Uses binary outputs from Step 2.1 and Step 2.2 as inputs. The significance threshold for et_sig is based on the number of variants after LD correction, found in the .log file for Step 1.__. The significance threshold for gw_sig is based on the total number of variants, which can be found in the .log file for Step 1.__. 
It is best to run this R script (and any R script) on the MSI OpenOnDemand platform (https://ondemand.msi.umn.edu/pun/sys/dashboard/) 
```
manhattan_plots.R
```






# VARIANT ANALYSIS BY CANDIDATE GENES

### Step 3.1 Subset general population VCF by cases and remove cases from general population VCF
```
sbatch subset_cases&genpop.slurm
```


### Step 3.2 Subset general population (excluding cases) and cases by candidate gene regions
```
sbatch subset_by_CG.slurm
```


### Step 3.3 Obtain the longest transcript for each variant
