# #Load Packages and Data ####
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("tidyverse")
# #install.packages("stringr")
# install.packages("readr")
# install.packages("VennDiagram")
# install.packages("devtools")
# devtools::install_github("gaospecial/ggVennDiagram")


# library(dplyr)
# library(tidyr)
# library(tidyverse)
# #library(stringr)
# library(readr)
# library(VennDiagram)
# 
# 
# library(BiocManager)
# # BiocManager::install("biomaRt")
# library(biomaRt)

pheno <- read_delim("/home/durwa004/stock526/CandidateGenes/data/Phenolyzer_june14_SeedGeneList.txt")
clinv <- read_delim("/home/durwa004/stock526/CandidateGenes/data/clinvar_result_june15_.txt")
OMIM_ca <- read_delim("/home/durwa004/stock526/CandidateGenes/data/OMIM_cardiac_arrhythmia_july26.tsv", skip = 4)
OMIM_sad <- read_delim("/home/durwa004/stock526/CandidateGenes/data/OMIM_SAD_july26.tsv", skip = 4)
OMIM_scd <- read_delim("/home/durwa004/stock526/CandidateGenes/data/OMIM_sudden_cardiac_death_july26.tsv", skip = 4)
OT_ca <- read_delim("/home/durwa004/stock526/CandidateGenes/data/open_targets_cardiac_arrhythmia_july26.tsv")
OT_sca <- read_delim("/home/durwa004/stock526/CandidateGenes/data/open_targets_sudden_cardiac_arrest_july26.tsv")

# OMIM data wrangling####

process_OMIM <- function(data_table) {
  # Extract the third column
  column <- as.character(data_table[, 3])
  
  # Split the column values by punctuation characters
  col_split <- unlist(lapply(column, function(x) {
    strsplit(x, "[[:punct:]]")
  }))
  
  # Convert the split values to character type and remove leading/trailing whitespace
  col_split <- trimws(as.character(col_split))
  
  # Remove any blank rows
  col_split <- na.omit(col_split)
  
  # Convert the result to a data frame
  result_table <- data.frame(Value = col_split)
  
  # Remove blank rows from the data frame
  result_table <- result_table[result_table$Value != "", , drop = FALSE]
  result_table<- unique(result_table)
  return(result_table)
}

OMIM_ca <- process_OMIM(OMIM_ca)
OMIM_sad <- process_OMIM(OMIM_sad)
OMIM_scd <- process_OMIM(OMIM_scd)

OMIM<-unique(union(union(OMIM_ca,OMIM_sad),OMIM_scd))

# OT data wrangling ####

process_OT <- function(data1) {
  data1 <- unique(data1[,1])
  names(data1)[names(data1) == "symbol"] <- "Value"
  return(data1)
}

OT_ca<-process_OT(OT_ca)
OT_sca <- process_OT(OT_sca)
OT<-unique(union(OT_ca,OT_sca))


# Phenolyzer data wrangling ####
pheno <- unique(pheno[,2])
names(pheno)[names(pheno) == "Gene"] <- "Value"


# ClinVar data wrangling ####

process_clinv <- function(data_table) {
  column <- as.character(data_table[, 2])
  
  # Split the column values by both ", " and "|"
  col_split <- unlist(lapply(column, function(x) {
    strsplit(x, ", |\\|")
  }))
  
  # Remove any blank rows
  col_split <- na.omit(col_split)
  
  # Remove all punctuation except hyphens
  col_split <- gsub("[\"\\(\\)]", "", col_split, perl = TRUE)
  
  # Convert the split values to character type and remove leading/trailing whitespace
  col_split <- trimws(as.character(col_split))
  
  # Convert the result to a data frame
  result_table <- data.frame(Value = col_split)
  
  # Remove blank rows from the data frame
  result_table <- result_table[result_table$Value != "", , drop = FALSE]
  
  result_table<-unique(result_table)
  
  return(result_table)
}


clinv <- process_clinv(clinv)

# Create Venn Diagrams ####


venn.diagram(
  x = c(OT = OT, OMIM = OMIM, ClinVar = clinv,  pheno = pheno),
  category.names = c("OT","OMIM", "ClinVar","Phenolyzer" ),
  filename = "four_set_venn.png",
  output = TRUE,
  imagetype = "png",
  scaled = TRUE,
  col = "black",
  fill = c("blue", "green", "red", "purple"),
  alpha = 0.5,
  label.col = "black",
  cex = 2,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.col = c("blue", "green", "red", "purple"),
  cat.cex = 2,
  cat.pos = c(0.2,0.2, 0.2, 0.2),  # Adjust cat.pos for all categories
  cat.dist = c(0.2, 0.2, 0.1, 0.1),  # Adjust cat.dist for all categories
  margin = 0.05,  # Increase overall margin
  spacing = 0.25  # Adjust spacing between circles
)
# Intersection ####


# Find narrow intersection
narrow <- unique(intersect(OMIM,intersect(OT,intersect(pheno,clinv))))


# Union of programs with multiple vectors
A <- clinv
B <- OMIM
C <- OT
D <- pheno

AB <- (unique(intersect( A,  B)))
AC <- (unique(intersect( A,  C)))
AD <- (unique(intersect( A,  D)))
BC <- (unique(intersect( B,  C)))
BD <- (unique(intersect( B,  D)))
CD <- (unique(intersect( C,  D)))

ACD <- (unique(intersect( A,  CD)))
BCD <- (unique(intersect( B,  CD)))
ABD <- (unique(intersect( A,  BD)))
ABC <- (unique(intersect( A,  BC)))

ABCD <- unique(intersect(OMIM,intersect(OT,intersect(pheno,clinv))))
  
broad <- Reduce(union, list(AB, AC, AD, BC, BD, CD, ACD, BCD, ABD, ABC))
broad <- unique(broad)  #692 genes




# BioMart Annotations ####

ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "ecaballus_gene_ensembl",
                      mirror = "useast")

# Retrieve all Genes in EqCab

BioMartData<- getBM( attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                     mart = ensembl)

#Rename BioMartData column gene name so that it matches the CG lists
names(BioMartData)[1] <- "Value"



# find data that has NAs when converting from human candidate genes to the horse genes
# temp<-left_join(input, BioMartData, by = "Value")
# I tried to make code for this but my brain was fried, just look at the data and write it down


# Manual annotations ####


# Manually annotate the orthologs so that they can be merged back together

# Because CALM2 is one of the narrow candidate genes (mined in the human genome),
# yet it doesnt have a proper ortholog in the horse genome I am replacing the data
# for CALM2 for CALM1 and CALM3 so that a narrow candidate gene is not lost. Additionally
# these two CALM genes will be added back into broad and narrow as well to make things 
# even. 

## scratch that, do it for all candidate genes that did not map to EquCCab 



# CALM1
c1 <- "CALM1"
CALM1 <- BioMartData[BioMartData$Value == c1, ]
c3 <- "CALM3"
CALM3 <- BioMartData[BioMartData$Value == c3, ]
CALM <- full_join(CALM1,CALM3)

# PHYH
# equine ortholog of PHYH is ENSECAG00000029564
ph <- "ENSECAG00000029564"
PHYH <- BioMartData[BioMartData$ensembl_gene_id == ph, ]
PHYH[1,1] <- "ENSECAG00000029564"

# SCN11A
# Equine ortholog of SCN11SA = GORASP1
# (https://useast.ensembl.org/Equus_caballus/Gene/Summary?g=ENSECAG00000016287;r=16:47909032-48054877)
sc <- "GORASP1"
SCN11A <- BioMartData[BioMartData$Value == sc, ]

# BAG3
# Does not have an equine ortholog

# TTN-AS1
# No equine ortholog, antisense RNA codeing 
# TTN is already included in the gene lists

# HAND2-AS1
# No equine ortholog, antisense RNA
# HAND2 is already in gene lists


#Inner_join candidate gene files and biomart data

# 1 --> has ensembl gene IDs
# 2 --> does not have ensembl gene IDs and start/stop are +- 5000 bp

# Use this function to annotate the Ensembl Gene IDs, start/stop positions,
# and lastly all of the genes present in the CG list but do not have a corresponding 
# ortholog in horses (eg, CALM1)



BioMart_annotations_1 <- function(input) {
  # Perform left join with BioMartData
  output1 <- left_join(input, BioMartData, by = "Value")
  
  # Full join
  output1 <- full_join(CALM, output1, by = c("Value", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"))
  output1 <- full_join(PHYH, output1, by = c("Value", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"))
  output1 <- full_join(SCN11A, output1, by = c("Value", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"))
  
  # Remove rows with NAs
  output1 <- output1[!is.na(output1[, 2]), ]
  
  # Exclude 'ensembl_gene_id' and 'Value' columns from broad1
  output2 <- output1[, !names(output1) %in% c("ensembl_gene_id", "Value"), drop = FALSE]
  
  # Append "chr" to the values in the first column of output2
  output2[, 1] <- paste0("chr", output2[, 1])
  
  # Add 5000 to start_position and subtract 5000 from end_position
  output2$start_position <- as.integer(output2$start_position - 5000)
  output2$end_position <- as.integer(output2$end_position + 5000)
  
  # Return gene1 and gene2
  return(list(gene1 = output1, gene2 = output2))
}


# broad candidate genes
names(broad)[1] <- "Value"
result <- BioMart_annotations_1(broad)
broad1 <- (result$gene1)
broad2 <- (result$gene2)

# narrow candidate genes
result <- BioMart_annotations_1(narrow)
narrow1 <- result$gene1
narrow2 <- result$gene2


# Export Final Results ####

# 1 --> has ensembl gene IDs
# 2 --> does not have ensembl gene IDs and start/stop are +- 5000 bp

save_vector <- function(my_vector, file_path) {
  # Write vector to a text file with tab separator and without header
  write_tsv(my_vector, file = file_path, col_names = FALSE)
}

save_vector(broad2, "/home/durwa004/stock526/CandidateGenes/results/broad_candidate_genes.txt")
save_vector(narrow2, "/home/durwa004/stock526/CandidateGenes/results/narrow_candidate_genes.txt")

