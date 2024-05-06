### Filter SNPs for doing the GWAS

##Running this on the CSC cluster
# Setting the paths
.libPaths(c("/projappl/project_2000350/rpackages", .libPaths()))
libpath <- .libPaths()[1]

#Load sources and libraries
source("/projappl/project_2000350/scripts/Neuro_functions.R") #This works
source("/projappl/project_2000350/gwas/association_scripts.R")
library(dplyr)

#Load the genotypes 
myG <- read.table(file = "/projappl/project_2000350/gwas/data/Neuro_hapmap_natpop_all_filtered.txt", header = FALSE, stringsAsFactors = FALSE)
#geno.names <- myG[1,-c(1:11)]

#This file contains strains from family E that were not included in this study, and have to be removed
#V376 - V426

myG <- myG[,-c(376:426)] #Drop family E

### Then filter for low minor allele counts

#Need to filter again for SNP with low minor allele counts
ma.counts <- calc.alleles(as.matrix(myG[,12:ncol(myG)])) #Calculate minor allele counts (better than freq since totals are not the same for all)
check <- ma.counts > 5
check[1] <- T ## Set the first row as true so that genotype names etc are not removed
myG <- myG[check,] #Filter SNPs with low minor allele counts

### Save new gnotype file

write.table(myG, file = "/projappl/project_2000350/gwas/data/Neuro_hapmap_natpop_gwas_temp.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
