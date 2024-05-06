#Trying some analysis of .vcf files (which are small) with R

library(tidyr)
library(dplyr)
library(stringr) #For counting matches in strings
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

#library(doMC)
#library(foreach)

#Load scripts for dealing with RAD-seq data and .vcf files
source("~/Documents/tutkijatohtori/scripts/vcfscripts.R")

#Note for different families: The assumption is that code blocks are run separately for each family

### * Process the .vcf file (from GATK)

### ** Family A

#Read .vcf file into R
snpdataUG <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/famA/famA.vcf", comment.char = "", sep = "\t", header = T, skip = 52, stringsAsFactors = FALSE)

#head(separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")[,1:30])

#Split the info column into number of samples with data and allele frequency
snpdataUG <- separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")
colnames(snpdataUG)[8] <- "NS"
snpdataUG$NS <- apply(snpdataUG[,-c(1:21)], MARGIN = 1, NS.count)
snpdataUG$AN <- as.integer(gsub("AN=", "", snpdataUG$AN))
snpdataUG$DP <- as.integer(gsub("DP=", "", snpdataUG$DP))
temp1 <- gsub("AF=", "", snpdataUG$AF)
temp1 <- strsplit(temp1, ",")
temp1 <- as.numeric(sapply(temp1, '[[', 1))
snpdataUG$AF <- calc.MAF.UG(temp1) ##Calculates minor allele frequency

#Get rid of loci that have less than 20 individuals genotyped
snpdataUG <- filter(snpdataUG, NS >= 20)

#Extract genotypes
genotypesUG <- process.vcf.genotypesUG(snpdataUG, cols = 22:117)

### ** Family B

#Read .vcf file into R
snpdataUG <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/famB/famB.vcf", comment.char = "", sep = "\t", header = T, skip = 52, stringsAsFactors = FALSE)

#Split the info column into number of samples with data and allele frequency
snpdataUG <- separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")
colnames(snpdataUG)[8] <- "NS"
snpdataUG$NS <- apply(snpdataUG[,-c(1:21)], MARGIN = 1, NS.count)
snpdataUG$AN <- as.integer(gsub("AN=", "", snpdataUG$AN))
snpdataUG$DP <- as.integer(gsub("DP=", "", snpdataUG$DP))
temp1 <- gsub("AF=", "", snpdataUG$AF)
temp1 <- strsplit(temp1, ",")
temp1 <- as.numeric(sapply(temp1, '[[', 1))
snpdataUG$AF <- calc.MAF.UG(temp1) ##Calculates minor allele frequency

#Get rid of loci that have less than 10 individuals genotyped
snpdataUG <- filter(snpdataUG, NS >= 10)

#Extract genotypes
genotypesUG <- process.vcf.genotypesUG(snpdataUG, cols = 22:73)

### ** Family C

#Read .vcf file into R
snpdataUG <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/famC/famC.vcf", comment.char = "", sep = "\t", header = T, skip = 52, stringsAsFactors = FALSE)

#head(separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")[,1:30])

#Split the info column into number of samples with data and allele frequency
snpdataUG <- separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")
colnames(snpdataUG)[8] <- "NS"
snpdataUG$NS <- apply(snpdataUG[,-c(1:21)], MARGIN = 1, NS.count)
snpdataUG$AN <- as.integer(gsub("AN=", "", snpdataUG$AN))
snpdataUG$DP <- as.integer(gsub("DP=", "", snpdataUG$DP))
temp1 <- gsub("AF=", "", snpdataUG$AF)
temp1 <- strsplit(temp1, ",")
temp1 <- as.numeric(sapply(temp1, '[[', 1))
snpdataUG$AF <- calc.MAF.UG(temp1) ##Calculates minor allele frequency

#Get rid of loci that have less than 10 individuals genotyped
snpdataUG <- filter(snpdataUG, NS >= 10)

#Extract genotypes
genotypesUG <- process.vcf.genotypesUG(snpdataUG, cols = 22:73)

### ** Family D

#Read .vcf file into R
snpdataUG <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/famD/famD.vcf", comment.char = "", sep = "\t", header = T, skip = 52, stringsAsFactors = FALSE)

#head(separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")[,1:30])

#Split the info column into number of samples with data and allele frequency
snpdataUG <- separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")
colnames(snpdataUG)[8] <- "NS"
snpdataUG$NS <- apply(snpdataUG[,-c(1:21)], MARGIN = 1, NS.count)
snpdataUG$AN <- as.integer(gsub("AN=", "", snpdataUG$AN))
snpdataUG$DP <- as.integer(gsub("DP=", "", snpdataUG$DP))
temp1 <- gsub("AF=", "", snpdataUG$AF)
temp1 <- strsplit(temp1, ",")
temp1 <- as.numeric(sapply(temp1, '[[', 1))
snpdataUG$AF <- calc.MAF.UG(temp1) ##Calculates minor allele frequency

#Get rid of loci that have less than 10 individuals genotyped
snpdataUG <- filter(snpdataUG, NS >= 10)

#Extract genotypes
genotypesUG <- process.vcf.genotypesUG(snpdataUG, cols = 22:75)


### ** Family E

#Read .vcf file into R
snpdataUG <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/famE/famE.vcf", comment.char = "", sep = "\t", header = T, skip = 52, stringsAsFactors = FALSE)

#head(separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")[,1:30])

#Split the info column into number of samples with data and allele frequency
snpdataUG <- separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")
colnames(snpdataUG)[8] <- "NS"
snpdataUG$NS <- apply(snpdataUG[,-c(1:21)], MARGIN = 1, NS.count)
snpdataUG$AN <- as.integer(gsub("AN=", "", snpdataUG$AN))
snpdataUG$DP <- as.integer(gsub("DP=", "", snpdataUG$DP))
temp1 <- gsub("AF=", "", snpdataUG$AF)
temp1 <- strsplit(temp1, ",")
temp1 <- as.numeric(sapply(temp1, '[[', 1))
snpdataUG$AF <- calc.MAF.UG(temp1) ##Calculates minor allele frequency

#Get rid of loci that have less than 10 individuals genotyped
snpdataUG <- filter(snpdataUG, NS >= 10)

#Extract genotypes
genotypesUG <- process.vcf.genotypesUG(snpdataUG, cols = 22:74)

### ** Family G

#Read .vcf file into R
snpdataUG <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/famG/famG.vcf", comment.char = "", sep = "\t", header = T, skip = 52, stringsAsFactors = FALSE)

#head(separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")[,1:30])

#Split the info column into number of samples with data and allele frequency
snpdataUG <- separate(snpdataUG, col = "INFO", into = c("AC", "AF", "AN", "DP", "Del", "FS", "HS", "ML1", "ML2", "MQ", "MQ0", "QD", "SO"), sep = ";")
colnames(snpdataUG)[8] <- "NS"
snpdataUG$NS <- apply(snpdataUG[,-c(1:21)], MARGIN = 1, NS.count)
snpdataUG$AN <- as.integer(gsub("AN=", "", snpdataUG$AN))
snpdataUG$DP <- as.integer(gsub("DP=", "", snpdataUG$DP))
temp1 <- gsub("AF=", "", snpdataUG$AF)
temp1 <- strsplit(temp1, ",")
temp1 <- as.numeric(sapply(temp1, '[[', 1))
snpdataUG$AF <- calc.MAF.UG(temp1) ##Calculates minor allele frequency

#Get rid of loci that have less than 10 individuals genotyped
snpdataUG <- filter(snpdataUG, NS >= 10)

#Extract genotypes
genotypesUG <- process.vcf.genotypesUG(snpdataUG, cols = 22:92)


### * Filter UG SNPs that are bad

### ** Family A

koe2 <- avg.quality.depth.vcf.UG(genotypesUG)
koe3 <- hets.vcf.UG(genotypesUG)
koe4 <- allelic.depth.balance.UG(genotypesUG)

gq.filter <- koe2[,1] >= 30 #Keep only SNPs that have average genotype quality >= 30
dp.filter <- koe2[,2] > 10 #Keep only SNPs with > 10 mean reads across all inds
af.filter <- snpdataUG$AF > 0.2 #Keep only SNPs that have minor allele frequency > 0.2
adrat.filter <- log10(koe4[,4]) > -1 & log10(koe4[,4]) < 1 #Keep only SNPs that have balanced allelic depth ratio

#Checking that both parents are not missing (another can be inferred) or are the same
parent.filter <- genotypesUG$P10948[,1] != genotypesUG$P10886[,1]

#One (or both) parent is heterozygous
parenthet.filter <- is.het.UG(genotypesUG$P10948[,1]) == F & is.het.UG(genotypesUG$P10886[,1]) == F

#HET filter
het.filter <- koe3 < 0.1 #Proportion of heterozygotes in genotyped samples should be < 0.1

##Perform filtering
#Make final filter
filters <- cbind(gq.filter, dp.filter, adrat.filter, parent.filter, parenthet.filter, het.filter, af.filter) #, adrat.filter, ad.filter)
final.filter <- apply(filters, 1, all)

snpdataUG.filt <- snpdataUG[final.filter,] #Final SNPs (at the moment 14476)

##Extract genotypes
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:117)

###Additional genotype quality control

#Set all genotype calls with seq depth <= 5 to missing data
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6)
#Check do loci with multiple segregating alleles exist
alleles.res <- check.allele.numbers(genotypesUG.filt)
#summary(alleles.res) #No monomorphic or triallelic loci exists
indels <- is.indel(snpdataUG.filt) #indels, 1599 in total
snpdataUG.filt$ID[indels] <- "indel"
snpdataUG.filt$ID[!indels] <- "SNP" #Store information about whether locus is a SNP or indel

duplicates <- duplicated(paste(snpdataUG.filt[,1], snpdataUG.filt[,2], sep = "_")) #Check for duplicated genome locations
#snpdataUG.filt[duplicates,1:21] ##18 duplicated locations in total
##These are all locations where there is both a SNP and an indel (complex mutation?)
##Keeping only indels (as these are probably result of complex mutations)
dup.ind <- snpdataUG.filt$POS %in% snpdataUG.filt[duplicates,1:21]$POS
#snpdataUG.filt[dup.ind,1:21]
snp.ind <- snpdataUG.filt$ID == "SNP"
##Cases that are both duplicated and SNP (of the duplicated pair)
dupfiltind <- apply(data.frame(dup.ind,snp.ind), MARGIN = 1, all)


##Remove one pair of the duplicates and extracting genotypes again
snpdataUG.filt <- snpdataUG.filt[!dupfiltind,] #14458 loci left
#Also filtering again for loci where both parents are missing (because missing data was inserted at bad genotypes, this causes a few cases of both parents being missing)
parent.filter2 <- genotypesUG.filt$P10948[,1] != genotypesUG.filt$P10886[,1]
snpdataUG.filt <- snpdataUG.filt[parent.filter2,] #Filters out 13 loci, 14445 left
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:117)
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6) #Run this again


##Some plotting to check stuff
ggplot(snpdataUG.filt, aes(x = NS, y = AF)) +
    geom_point()

ggplot(snpdataUG, aes(x = NS, y = AF)) +
    geom_point()

ggplot(data.frame(koe4), aes(x = log10(ratio))) +
    geom_histogram()

### ** Family B

koe2 <- avg.quality.depth.vcf.UG(genotypesUG)
koe3 <- hets.vcf.UG(genotypesUG)
koe4 <- allelic.depth.balance.UG(genotypesUG)

gq.filter <- koe2[,1] >= 30 #Keep only SNPs that have average genotype quality >= 30
dp.filter <- koe2[,2] > 10 #Keep only SNPs with > 10 mean reads across all inds
af.filter <- snpdataUG$AF > 0.2 #Keep only SNPs that have minor allele frequency > 0.2
adrat.filter <- log10(koe4[,4]) > -1 & log10(koe4[,4]) < 1 #Keep only SNPs that have balanced allelic depth ratio

#Checking that both parents are not missing (another can be inferred) or are the same
parent.filter <- genotypesUG$P10932[,1] != genotypesUG$P1165[,1]

#One (or both) parent is heterozygous
parenthet.filter <- is.het.UG(genotypesUG$P10932[,1]) == F & is.het.UG(genotypesUG$P1165[,1]) == F

#HET filter
het.filter <- koe3 < 0.1 #Proportion of heterozygotes in genotyped samples should be < 0.1

##Perform filtering
#Make final filter
filters <- cbind(gq.filter, dp.filter, adrat.filter, parent.filter, parenthet.filter, het.filter, af.filter) #, adrat.filter, ad.filter)
final.filter <- apply(filters, 1, all)

snpdataUG.filt <- snpdataUG[final.filter,] #Number of SNPs before filtering 47486
#After filtering 15728

##Extract genotypes
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:73)

###Additional genotype quality control

#Set all genotype calls with seq depth <= 5 to missing data
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6)
#Check do loci with multiple segregating alleles exist
alleles.res <- check.allele.numbers(genotypesUG.filt)
#summary(alleles.res) #No monomorphic or triallelic loci exists
indels <- is.indel(snpdataUG.filt) #indels, 1769 in total
snpdataUG.filt$ID[indels] <- "indel"
snpdataUG.filt$ID[!indels] <- "SNP" #Store information about whether locus is a SNP or indel

duplicates <- duplicated(paste(snpdataUG.filt[,1], snpdataUG.filt[,2], sep = "_")) #Check for duplicated genome locations
#snpdataUG.filt[duplicates,1:21] ##28 duplicated locations in total
##These are all locations where there is both a SNP and an indel (complex mutation?)
##Keeping only indels (as these are probably result of complex mutations)
dup.ind <- snpdataUG.filt$POS %in% snpdataUG.filt[duplicates,1:21]$POS
#snpdataUG.filt[dup.ind,1:21]
snp.ind <- snpdataUG.filt$ID == "SNP"
##Cases that are both duplicated and SNP (of the duplicated pair)
dupfiltind <- apply(data.frame(dup.ind,snp.ind), MARGIN = 1, all)

#Also filtering again for loci where both parents are missing (because missing data was inserted at bad genotypes, this causes a few cases of both parents being missing)
parent.filter2 <- genotypesUG.filt$P10932[,1] != genotypesUG.filt$P1165[,1]

#Combining filters
filter2 <- cbind(!dupfiltind, parent.filter2)
final.filter2 <- apply(filter2, 1, all)

##Remove one pair of the duplicates and extracting genotypes again
snpdataUG.filt <- snpdataUG.filt[final.filter2,] #13523 loci left
#snpdataUG.filt <- snpdataUG.filt[-13616,] #Removed one bad locus, 15680 loci left
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:73)
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6) #Run this again

##Some plotting to check stuff
ggplot(snpdataUG.filt, aes(x = NS, y = AF)) +
    geom_point()

ggplot(snpdataUG, aes(x = NS, y = AF)) +
    geom_point()

ggplot(data.frame(koe4), aes(x = log10(ratio))) +
    geom_histogram()


### ** Family C

koe2 <- avg.quality.depth.vcf.UG(genotypesUG)
koe3 <- hets.vcf.UG(genotypesUG)
koe4 <- allelic.depth.balance.UG(genotypesUG)

gq.filter <- koe2[,1] >= 30 #Keep only SNPs that have average genotype quality >= 30
dp.filter <- koe2[,2] > 10 #Keep only SNPs with > 10 mean reads across all inds
af.filter <- snpdataUG$AF > 0.2 #Keep only SNPs that have minor allele frequency > 0.2
adrat.filter <- log10(koe4[,4]) > -1 & log10(koe4[,4]) < 1 #Keep only SNPs that have balanced allelic depth ratio

#Checking that both parents are not missing (another can be inferred) or are the same
parent.filter <- genotypesUG$P4498[,1] != genotypesUG$P8816[,1]

#One (or both) parent is heterozygous
parenthet.filter <- is.het.UG(genotypesUG$P4498[,1]) == F & is.het.UG(genotypesUG$P8816[,1]) == F

#HET filter
het.filter <- koe3 < 0.1 #Proportion of heterozygotes in genotyped samples should be < 0.1

##Perform filtering
#Make final filter
filters <- cbind(gq.filter, dp.filter, adrat.filter, parent.filter, parenthet.filter, het.filter, af.filter) #, adrat.filter, ad.filter)
final.filter <- apply(filters, 1, all)

snpdataUG.filt <- snpdataUG[final.filter,] #Number of SNPs before filtering 48198
#After filtering 17442

##Extract genotypes
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:73)

###Additional genotype quality control

#Set all genotype calls with seq depth <= 5 to missing data
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6)
#Check do loci with multiple segregating alleles exist
alleles.res <- check.allele.numbers(genotypesUG.filt)
#summary(alleles.res) #No monomorphic or triallelic loci exists
indels <- is.indel(snpdataUG.filt) #indels, 1769 in total
snpdataUG.filt$ID[indels] <- "indel"
snpdataUG.filt$ID[!indels] <- "SNP" #Store information about whether locus is a SNP or indel

duplicates <- duplicated(paste(snpdataUG.filt[,1], snpdataUG.filt[,2], sep = "_")) #Check for duplicated genome locations
#snpdataUG.filt[duplicates,1:21] ##17 duplicated locations in total
##These are all locations where there is both a SNP and an indel (complex mutation?)
##Keeping only indels (as these are probably result of complex mutations)
dup.ind <- snpdataUG.filt$POS %in% snpdataUG.filt[duplicates,1:21]$POS
#snpdataUG.filt[dup.ind,1:21]
snp.ind <- snpdataUG.filt$ID == "SNP"
##Cases that are both duplicated and SNP (of the duplicated pair)
dupfiltind <- apply(data.frame(dup.ind,snp.ind), MARGIN = 1, all)

#Also filtering again for loci where both parents are missing (because missing data was inserted at bad genotypes, this causes a few cases of both parents being missing)
parent.filter2 <- genotypesUG.filt$P4498[,1] != genotypesUG.filt$P8816[,1]

#Combining filters
filter2 <- cbind(!dupfiltind, parent.filter2)
final.filter2 <- apply(filter2, 1, all)

##Remove one pair of the duplicates and extracting genotypes again
snpdataUG.filt <- snpdataUG.filt[final.filter2,] #17297 loci left
#snpdataUG.filt <- snpdataUG.filt[-13616,] #Removed one bad locus, 15680 loci left
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:73)
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6) #Run this again

#Applying an extra allele frequency filter
#af.filter <- snpdataUG.filt$AF > 0.2
#snpdataUG.filt <- snpdataUG.filt[af.filter,]



##Some plotting to check stuff
ggplot(snpdataUG.filt, aes(x = NS, y = AF)) +
    geom_point()

ggplot(snpdataUG, aes(x = NS, y = AF)) +
    geom_point()

ggplot(data.frame(koe4), aes(x = log10(ratio))) +
    geom_histogram()


### ** Family D

koe2 <- avg.quality.depth.vcf.UG(genotypesUG)
koe3 <- hets.vcf.UG(genotypesUG)
koe4 <- allelic.depth.balance.UG(genotypesUG)

gq.filter <- koe2[,1] >= 30 #Keep only SNPs that have average genotype quality >= 30
dp.filter <- koe2[,2] > 10 #Keep only SNPs with > 10 mean reads across all inds
af.filter <- snpdataUG$AF > 0.2 #Keep only SNPs that have minor allele frequency > 0.2
adrat.filter <- log10(koe4[,4]) > -1 & log10(koe4[,4]) < 1 #Keep only SNPs that have balanced allelic depth ratio

#Checking that both parents are not missing (another can be inferred) or are the same
parent.filter <- genotypesUG$P3223[,1] != genotypesUG$P8845[,1]

#One (or both) parent is heterozygous
parenthet.filter <- is.het.UG(genotypesUG$P3223[,1]) == F & is.het.UG(genotypesUG$P8845[,1]) == F

#HET filter
het.filter <- koe3 < 0.1 #Proportion of heterozygotes in genotyped samples should be < 0.1

##Perform filtering
#Make final filter
filters <- cbind(gq.filter, dp.filter, adrat.filter, parent.filter, parenthet.filter, het.filter, af.filter) #, adrat.filter, ad.filter)
final.filter <- apply(filters, 1, all)

snpdataUG.filt <- snpdataUG[final.filter,] #Number of SNPs before filtering 43165
#After filtering 16352

##Extract genotypes
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:75)

###Additional genotype quality control

#Set all genotype calls with seq depth <= 5 to missing data
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6)
#Check do loci with multiple segregating alleles exist
alleles.res <- check.allele.numbers(genotypesUG.filt)
#summary(alleles.res) #No monomorphic or triallelic loci exists
indels <- is.indel(snpdataUG.filt) #indels, 1886 in total
snpdataUG.filt$ID[indels] <- "indel"
snpdataUG.filt$ID[!indels] <- "SNP" #Store information about whether locus is a SNP or indel

duplicates <- duplicated(paste(snpdataUG.filt[,1], snpdataUG.filt[,2], sep = "_")) #Check for duplicated genome locations
#snpdataUG.filt[duplicates,1:21] ##36 duplicated locations in total
##These are all locations where there is both a SNP and an indel (complex mutation?)
##Keeping only indels (as these are probably result of complex mutations)
dup.ind <- snpdataUG.filt$POS %in% snpdataUG.filt[duplicates,1:21]$POS
#snpdataUG.filt[dup.ind,1:21]
snp.ind <- snpdataUG.filt$ID == "SNP"
##Cases that are both duplicated and SNP (of the duplicated pair)
dupfiltind <- apply(data.frame(dup.ind,snp.ind), MARGIN = 1, all)

#Also filtering again for loci where both parents are missing (because missing data was inserted at bad genotypes, this causes a few cases of both parents being missing)
parent.filter2 <- genotypesUG.filt$P3223[,1] != genotypesUG.filt$P8845[,1]

#Combining filters
filter2 <- cbind(!dupfiltind, parent.filter2)
final.filter2 <- apply(filter2, 1, all)

##Remove one pair of the duplicates and extracting genotypes again
snpdataUG.filt <- snpdataUG.filt[final.filter2,] #16314 loci left
snpdataUG.filt <- snpdataUG.filt[-(16312:16314),] #Removed three mtDNA markers, 16311 loci left

genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:75)
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6) #Run this again

##Some plotting to check stuff
ggplot(snpdataUG.filt, aes(x = NS, y = AF)) +
    geom_point()

ggplot(snpdataUG, aes(x = NS, y = AF)) +
    geom_point()

ggplot(data.frame(koe4), aes(x = log10(ratio))) +
    geom_histogram()



### ** Family E

koe2 <- avg.quality.depth.vcf.UG(genotypesUG)
koe3 <- hets.vcf.UG(genotypesUG)
koe4 <- allelic.depth.balance.UG(genotypesUG)

gq.filter <- koe2[,1] >= 30 #Keep only SNPs that have average genotype quality >= 30
dp.filter <- koe2[,2] > 10 #Keep only SNPs with > 10 mean reads across all inds
af.filter <- snpdataUG$AF > 0.2 #Keep only SNPs that have minor allele frequency > 0.2
adrat.filter <- log10(koe4[,4]) > -1 & log10(koe4[,4]) < 1 #Keep only SNPs that have balanced allelic depth ratio

#Checking that both parents are not missing (another can be inferred) or are the same
parent.filter <- genotypesUG$P10908[,1] != genotypesUG$P847[,1]

#One (or both) parent is heterozygous
parenthet.filter <- is.het.UG(genotypesUG$P10908[,1]) == F & is.het.UG(genotypesUG$P847[,1]) == F

#HET filter
het.filter <- koe3 < 0.1 #Proportion of heterozygotes in genotyped samples should be < 0.1

##Perform filtering
#Make final filter
filters <- cbind(gq.filter, dp.filter, adrat.filter, parent.filter, parenthet.filter, het.filter, af.filter) #, adrat.filter, ad.filter)
final.filter <- apply(filters, 1, all)

snpdataUG.filt <- snpdataUG[final.filter,] #Number of SNPs before filtering 41327
#After filtering 13712

##Extract genotypes
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:74)

###Additional genotype quality control

#Set all genotype calls with seq depth <= 5 to missing data
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6)
#Check do loci with multiple segregating alleles exist
alleles.res <- check.allele.numbers(genotypesUG.filt)
#summary(alleles.res) #No monomorphic or triallelic loci exists
indels <- is.indel(snpdataUG.filt) #indels, 1434 in total
snpdataUG.filt$ID[indels] <- "indel"
snpdataUG.filt$ID[!indels] <- "SNP" #Store information about whether locus is a SNP or indel

duplicates <- duplicated(paste(snpdataUG.filt[,1], snpdataUG.filt[,2], sep = "_")) #Check for duplicated genome locations
#snpdataUG.filt[duplicates,1:21] ##13 duplicated locations in total
##These are all locations where there is both a SNP and an indel (complex mutation?)
##Keeping only indels (as these are probably result of complex mutations)
dup.ind <- snpdataUG.filt$POS %in% snpdataUG.filt[duplicates,1:21]$POS
#snpdataUG.filt[dup.ind,1:21]
snp.ind <- snpdataUG.filt$ID == "SNP"
##Cases that are both duplicated and SNP (of the duplicated pair)
dupfiltind <- apply(data.frame(dup.ind,snp.ind), MARGIN = 1, all)

#Also filtering again for loci where both parents are missing (because missing data was inserted at bad genotypes, this causes a few cases of both parents being missing)
parent.filter2 <- genotypesUG.filt$P10908[,1] != genotypesUG.filt$P847[,1]

#Combining filters
filter2 <- cbind(!dupfiltind, parent.filter2)
final.filter2 <- apply(filter2, 1, all)

##Remove one pair of the duplicates and extracting genotypes again
snpdataUG.filt <- snpdataUG.filt[final.filter2,] #13697 loci left
#Need to deal with one genotyping error that causes problems
snpdataUG.filt[7075,59] <- "./.:0,0,6:6:18:253,253,253,18,18,0"

genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:74)
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6) #Run this again

##Some plotting to check stuff
ggplot(snpdataUG.filt, aes(x = NS, y = AF)) +
    geom_point()

ggplot(snpdataUG, aes(x = NS, y = AF)) +
    geom_point()

ggplot(data.frame(koe4), aes(x = log10(ratio))) +
    geom_histogram()



### ** Family G

koe2 <- avg.quality.depth.vcf.UG(genotypesUG)
koe3 <- hets.vcf.UG(genotypesUG)
koe4 <- allelic.depth.balance.UG(genotypesUG)

gq.filter <- koe2[,1] >= 30 #Keep only SNPs that have average genotype quality >= 30
dp.filter <- koe2[,2] > 10 #Keep only SNPs with > 10 mean reads across all inds
af.filter <- snpdataUG$AF > 0.2 #Keep only SNPs that have minor allele frequency > 0.2
adrat.filter <- log10(koe4[,4]) > -1 & log10(koe4[,4]) < 1 #Keep only SNPs that have balanced allelic depth ratio

#Checking that both parents are not missing (another can be inferred) or are the same
parent.filter <- genotypesUG$P10904[,1] != genotypesUG$P851[,1]

#One (or both) parent is heterozygous
parenthet.filter <- is.het.UG(genotypesUG$P10904[,1]) == F & is.het.UG(genotypesUG$P851[,1]) == F

#HET filter
het.filter <- koe3 < 0.1 #Proportion of heterozygotes in genotyped samples should be < 0.1

##Perform filtering
#Make final filter
filters <- cbind(gq.filter, dp.filter, adrat.filter, parent.filter, parenthet.filter, het.filter, af.filter) #, adrat.filter, ad.filter)
final.filter <- apply(filters, 1, all)

snpdataUG.filt <- snpdataUG[final.filter,] #Number of SNPs before filtering 54008
#After filtering 11284

##Extract genotypes
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:92)

###Additional genotype quality control

#Set all genotype calls with seq depth <= 5 to missing data
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6)
#Check do loci with multiple segregating alleles exist
alleles.res <- check.allele.numbers(genotypesUG.filt)
#summary(alleles.res) #3 monomorphic loci
monofilt <- !alleles.res[,2] #
snpdataUG.filt <- snpdataUG.filt[monofilt,] 
indels <- is.indel(snpdataUG.filt) #indels, 1209 in total
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:92) #Need to include this here
snpdataUG.filt$ID[indels] <- "indel"
snpdataUG.filt$ID[!indels] <- "SNP" #Store information about whether locus is a SNP or indel

duplicates <- duplicated(paste(snpdataUG.filt[,1], snpdataUG.filt[,2], sep = "_")) #Check for duplicated genome locations
#snpdataUG.filt[duplicates,1:21] ##20 duplicated locations in total
##These are all locations where there is both a SNP and an indel (complex mutation?)
##Keeping only indels (as these are probably result of complex mutations)
dup.ind <- snpdataUG.filt$POS %in% snpdataUG.filt[duplicates,1:21]$POS
#snpdataUG.filt[dup.ind,1:21]
snp.ind <- snpdataUG.filt$ID == "SNP"
##Cases that are both duplicated and SNP (of the duplicated pair)
dupfiltind <- apply(data.frame(dup.ind,snp.ind), MARGIN = 1, all)

#Also filtering again for loci where both parents are missing (because missing data was inserted at bad genotypes, this causes a few cases of both parents being missing)
parent.filter2 <- genotypesUG.filt$P10904[,1] != genotypesUG.filt$P851[,1]

#Combining filters
filter2 <- cbind(!dupfiltind, parent.filter2)
final.filter2 <- apply(filter2, 1, all)

##Remove one pair of the duplicates and extracting genotypes again
snpdataUG.filt <- snpdataUG.filt[final.filter2,] #11261 loci left
#Need to deal with one genotyping error that causes problems
#snpdataUG.filt[7075,59] <- "./.:0,0,6:6:18:253,253,253,18,18,0"

genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:92)
genotypesUG.filt <- quality.control.genot(genotypesUG.filt, min.seq.depth = 6) #Run this again

#Because quality control filtering was done and that removes some stuff need to run this again
#Checking that both parents are not missing (another can be inferred) or are the same
parent.filter <- genotypesUG.filt$P10904[,1] != genotypesUG.filt$P851[,1]
snpdataUG.filt <- snpdataUG.filt[parent.filter,] #11225 loci left
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:92)
alleles.res <- check.allele.numbers(genotypesUG.filt)
bial.filter <- !alleles.res[,3]
snpdataUG.filt <- snpdataUG.filt[bial.filter,] #11162 loci left
genotypesUG.filt <- process.vcf.genotypesUG(snpdataUG.filt, cols = 22:92)

##Some plotting to check stuff
ggplot(snpdataUG.filt, aes(x = NS, y = AF)) +
    geom_point()

ggplot(snpdataUG, aes(x = NS, y = AF)) +
    geom_point()

ggplot(data.frame(koe4), aes(x = log10(ratio))) +
    geom_histogram()

### * Genotype calls to SNPs and parent calls (UG format)

### ** Family A

#Make a SNP column for each individual
genotypesUG.filt <- vcf.gts2snps.UG(snpdataUG.filt, genotypesUG.filt)

##Make a genotype column indicating from which parent a site is from
##For data input into R/qtl etc.
genotypesUG.filt <- vcf.parentgenot(genotypesUG.filt, 95, 96)

#save(snpdataUG.filt, genotypesUG.filt, file = "~/Genomics/Neurospora/RAD/genotypes/famA/famAgenotypes.RData")

#Can load famA genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famA/famAgenotypes.RData")

##Exporting data into R/qtl | ASMap format
famA.asmap <- export.vcf2asmap(snpdataUG.filt, genotypesUG.filt)
write.table(famA.asmap, file = "~/Genomics/Neurospora/RAD/genotypes/famA/famA_asmap.txt", quote = FALSE, sep = "\t", row.names = TRUE)

##Plotting to check some stuff
#plot allele frequency along the chromosome
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_line() +
    facet_grid(X.CHROM ~ .)

##Plotting points shows better average local allele freq of markers
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_point() +
    ylab("Minor allele frequency") +
    scale_y_continuous(limits = c(0, 0.5)) +
    facet_grid(X.CHROM ~ .)

### ** Family B

#Make a SNP column for each individual
genotypesUG.filt <- vcf.gts2snps.UG(snpdataUG.filt, genotypesUG.filt)

##Make a genotype column indicating from which parent a site is from
##For data input into R/qtl etc.
genotypesUG.filt <- vcf.parentgenot(genotypesUG.filt, 51, 52)

#save(snpdataUG.filt, genotypesUG.filt, file = "~/Genomics/Neurospora/RAD/genotypes/famB/famBgenotypes.RData")

#Can load famA genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famB/famBgenotypes.RData")

##Exporting data into R/qtl | ASMap format
famB.asmap <- export.vcf2asmap(snpdataUG.filt, genotypesUG.filt)
write.table(famB.asmap, file = "~/Genomics/Neurospora/RAD/genotypes/famB/famB_asmap.txt", quote = FALSE, sep = "\t", row.names = TRUE)

##Plotting to check some stuff
#plot allele frequency along the chromosome
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_line() +
    facet_grid(X.CHROM ~ .)

##Plotting points shows better average local allele freq of markers
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_point() +
    ylab("Minor allele frequency") +
    scale_y_continuous(limits = c(0, 0.5)) +
    facet_grid(X.CHROM ~ .)

### ** Family C

#Make a SNP column for each individual
genotypesUG.filt <- vcf.gts2snps.UG(snpdataUG.filt, genotypesUG.filt)

##Make a genotype column indicating from which parent a site is from
##For data input into R/qtl etc.
genotypesUG.filt <- vcf.parentgenot(genotypesUG.filt, 51, 52)

#save(snpdataUG.filt, genotypesUG.filt, file = "~/Genomics/Neurospora/RAD/genotypes/famC/famCgenotypes.RData")

#Can load famC genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famC/famCgenotypes.RData")

##Exporting data into R/qtl | ASMap format
famC.asmap <- export.vcf2asmap(snpdataUG.filt, genotypesUG.filt)
write.table(famC.asmap, file = "~/Genomics/Neurospora/RAD/genotypes/famC/famC_asmap.txt", quote = FALSE, sep = "\t", row.names = TRUE)

##Plotting to check some stuff
#plot allele frequency along the chromosome
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_line() +
    facet_grid(X.CHROM ~ .)

##Plotting points shows better average local allele freq of markers
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_point() +
    ylab("Minor allele frequency") +
    scale_y_continuous(limits = c(0, 0.5)) +
    facet_grid(X.CHROM ~ .)



### ** Family D

#Make a SNP column for each individual
genotypesUG.filt <- vcf.gts2snps.UG(snpdataUG.filt, genotypesUG.filt)

##Make a genotype column indicating from which parent a site is from
##For data input into R/qtl etc.
genotypesUG.filt <- vcf.parentgenot(genotypesUG.filt, 53, 54)

#save(snpdataUG.filt, genotypesUG.filt, file = "~/Genomics/Neurospora/RAD/genotypes/famD/famDgenotypes.RData")

#Can load famA genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famD/famDgenotypes.RData")

##Exporting data into R/qtl | ASMap format
famD.asmap <- export.vcf2asmap(snpdataUG.filt, genotypesUG.filt)
write.table(famD.asmap, file = "~/Genomics/Neurospora/RAD/genotypes/famD/famD_asmap.txt", quote = FALSE, sep = "\t", row.names = TRUE)

##Plotting to check some stuff
#plot allele frequency along the chromosome
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_line() +
    facet_grid(X.CHROM ~ .)

##Plotting points shows better average local allele freq of markers
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_point() +
    ylab("Minor allele frequency") +
    scale_y_continuous(limits = c(0, 0.5)) +
    facet_grid(X.CHROM ~ .)

### ** Family E

#Make a SNP column for each individual
genotypesUG.filt <- vcf.gts2snps.UG(snpdataUG.filt, genotypesUG.filt)

##Make a genotype column indicating from which parent a site is from
##For data input into R/qtl etc.
genotypesUG.filt <- vcf.parentgenot(genotypesUG.filt, 52, 53)

#save(snpdataUG.filt, genotypesUG.filt, file = "~/Genomics/Neurospora/RAD/genotypes/famE/famEgenotypes.RData")

#Can load famA genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famE/famEgenotypes.RData")

##Exporting data into R/qtl | ASMap format
famE.asmap <- export.vcf2asmap(snpdataUG.filt, genotypesUG.filt)
write.table(famE.asmap, file = "~/Genomics/Neurospora/RAD/genotypes/famE/famE_asmap.txt", quote = FALSE, sep = "\t", row.names = TRUE)

##Plotting to check some stuff
#plot allele frequency along the chromosome
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_line() +
    facet_grid(X.CHROM ~ .)

##Plotting points shows better average local allele freq of markers
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_point() +
    ylab("Minor allele frequency") +
    scale_y_continuous(limits = c(0, 0.5)) +
    facet_grid(X.CHROM ~ .)

### ** Family G


#Make a SNP column for each individual
genotypesUG.filt <- vcf.gts2snps.UG(snpdataUG.filt, genotypesUG.filt)

##Make a genotype column indicating from which parent a site is from
##For data input into R/qtl etc.
genotypesUG.filt <- vcf.parentgenot(genotypesUG.filt, 70, 71)

#save(snpdataUG.filt, genotypesUG.filt, file = "~/Genomics/Neurospora/RAD/genotypes/famG/famGgenotypes.RData")

#Can load famA genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famG/famGgenotypes.RData")

##Exporting data into R/qtl | ASMap format
famG.asmap <- export.vcf2asmap(snpdataUG.filt, genotypesUG.filt)
write.table(famG.asmap, file = "~/Genomics/Neurospora/RAD/genotypes/famG/famG_asmap.txt", quote = FALSE, sep = "\t", row.names = TRUE)

##Plotting to check some stuff
#plot allele frequency along the chromosome
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_line() +
    facet_grid(X.CHROM ~ .)

##Plotting points shows better average local allele freq of markers
ggplot(snpdataUG.filt, aes(x = POS, y = AF)) +
    geom_point() +
    ylab("Minor allele frequency") +
    scale_y_continuous(limits = c(0, 0.5)) +
    facet_grid(X.CHROM ~ .)

### * Inferring genotypes for offspring based on parents

### ** Making the inferred combined SNP dataset

##Load all the SNPs from natural population reseuencing
snp1 <- read.table(file = "~/Genomics/Neurospora/natpop/Neuro_hapmap_natpop.txt", header = T, stringsAsFactors = FALSE)

##Drop all unmapped contigs
snp1 <- filter(snp1, chrom < 8)

##Load these one at a time and run the code below for each family
#Can load famA genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famA/famAgenotypes.RData")
#Can load famB genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famB/famBgenotypes.RData")
#Can load famC genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famC/famCgenotypes.RData")
#Can load famD genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famD/famDgenotypes.RData")
#Can load famE genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famE/famEgenotypes.RData")
#Can load famG genotypes
#load(file = "~/Genomics/Neurospora/RAD/genotypes/famG/famGgenotypes.RData")

##Run the code below for each family separately

##Modify snpdataUG.filt, so that chromosome is OK
snpdataUG.filt[,1] <- as.numeric(gsub("Supercontig_12.", "", snpdataUG.filt[,1]))
colnames(snpdataUG.filt)[1] <- "CHROM"

##Remember to change the parents between the families!
#Family A
#parentA <- "10886"
#parentB <- "10948"
#Family B
#parentA <- "10932"
#parentB <- "1165"
#Family C
#parentA <- "4498"
#parentB <- "8816"
#Family D
#parentA <- "3223"
#parentB <- "8845"
#Family E
#parentA <- "10908"
#parentB <- "847"
#Family G
parentA <- "10904"
parentB <- "851"

##Check that parent are OK
#head(genotypesUG.filt[['P10886']]) #X10886 is parent A
#head(genotypesUG.filt[['P10948']]) #X10948 is parent B

##Loop over all individuals (assumes parents are the last two are the parents that are dropped
fam.names <- names(genotypesUG.filt)[1:(length(names(genotypesUG.filt))-2)]

for(i in 1:length(fam.names)) {
    
    current.ind <- genotypesUG.filt[[i]] #Take current ind
    currentname <- names(genotypesUG.filt)[i] #Name of current inf
    
    SNPs <- NULL #Initialize variable
    
    ##Loop over all chromosomes
    for(j in 1:7) {
        sel.index <- snpdataUG.filt$CHROM == j
        current.chrom <- current.ind[sel.index,] #Genotypes, filtered by chrom

        segA <- find.segments(current.chrom$GENOP, "A")
        segB <- find.segments(current.chrom$GENOP, "B")

        genomic.coords <- segments.2.genomic.coordinates(segA, segB, snpdataUG.filt[sel.index,], chr = j)
        chrom.snps <- filter(snp1, chrom == j) #Take only current chromosome
        current.chrom <- infer.SNPs(genomicA = genomic.coords$A, genomicB = genomic.coords$B, parentA, parentB, currentname = currentname, allsnps = chrom.snps)
        SNPs <- rbind(SNPs, current.chrom) #Combine, should already be sorted by position
        
    }

    #Store the genotypes
    finalname <- paste("X", currentname, sep = "")
    snp1[[finalname]] <- SNPs[,2]

}


##Now snp1 contains the genotypes of natural strains and all families
##Save the gnotype file
write.table(snp1, file = "Neuro_hapmap_combined_all.txt", quote = FALSE, sep = "\t", row.names = FALSE)


### *** Making a plot of individual genotypes


plotdata <- data.frame(CHROM = snpdataUG.filt$CHROM, POS = snpdataUG.filt$POS, GENO = genotypesUG.filt$A2$GENOP)

pdf(file = "A2_genotype.pdf", height = 7, width = 7*1.618)
ggplot(plotdata, aes(x = POS, y = 1, color = GENO)) +
    geom_point() +
    #scale_x_continuous(expand = c(10,0)) +
    xlab("Position (bp)") +
    theme(axis.line.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none") +
    facet_grid(CHROM ~ .)
dev.off()    

### ** Filtering the combined datasset

library(doMC)
library(foreach)
registerDoMC(2)


#source("~/Documents/tutkijatohtori/association/association_scripts.R")

#Some functions that are needed
##Here use is.biallelic2(as.matrix(input)). This speeds things up considerably
is.biallelic2 <- function(aineisto) {
  results <- rep(0, dim(aineisto)[1])
  for(i in 1:dim(aineisto)[1]) {
      results[i] <- length(table(aineisto[i,], exclude = "N")) == 2
  }
  return(results)
}

#Calculate minor allele frequency, use: calc.MF(as.matrix(input))
calc.MAF <- function(aineisto) {
    results <- rep(0, dim(aineisto)[1])
    for(i in 1:dim(aineisto)[1]) {
        count <- table(aineisto[i,], exclude = "N")
        freqs <- c(count[1]/sum(count), count[2]/sum(count))
        results[i] <- min(freqs)
    }
    return(results)
}

#Assuming that snp1 contains all genotypes
snp1 <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/Neuro_hapmap_combined_all.txt", header = T, stringsAsFactors = FALSE)

### *** Performing some hacks, because not enough memory...

#Split this thing into two
write.table(snp1[1:500000,], file = "Neuro_combined_part1.txt", quote = FALSE, sep = "\t", row.names = FALSE)

write.table(snp1[500001:dim(snp1)[1],], file = "Neuro_combined_part2.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#### PART 1 ##############
snp1 <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/Neuro_combined_part1.txt", header = T, stringsAsFactors = FALSE)

##Check that SNPs are biallelic
testi <- is.biallelic2(as.matrix(snp1[,-c(1:11)])) #Note! Need to use as.matrix() here, otherwise too slow
#Now this running this function takes ~4 min

snp1 <- snp1[testi == 1,] ##531 SNPs that were not biallelic

##Calculating minor allele frequency
maf1 <- calc.MAF(as.matrix(snp1[,-c(1:11)]))

##Applying a minor allele frequency filter of 0.008
snp1 <- snp1[maf1 > 0.008,]

write.table(snp1, file = "Neuro_filtered_part1.txt", quote = FALSE, sep = "\t", row.names = FALSE) 
##############################

#### PART 2 #####################
snp2 <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/Neuro_combined_part2.txt", header = T, stringsAsFactors = FALSE)

testi <- is.biallelic2(as.matrix(snp2[,-c(1:11)]))

snp2 <- snp2[testi == 1,]

maf2 <- calc.MAF(as.matrix(snp2[,-c(1:11)]))

snp2 <- snp2[maf2 > 0.008,]

write.table(snp2, file = "Neuro_filtered_part2.txt", quote = FALSE, sep = "\t", row.names = FALSE) 
#################################

##### Combine the files ###############

part1 <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/Neuro_filtered_part1.txt", header = T, stringsAsFactors = FALSE)

part2 <- read.table(file = "~/Genomics/Neurospora/RAD/genotypes/Neuro_filtered_part2.txt", header = T, stringsAsFactors = FALSE)

snp <- rbind(part1, part2)

##Saving the combined filtered dataset
write.table(snp, file = "Neuro_hapmap_combined_all_filtered.txt", quote = FALSE, sep = "\t", row.names = FALSE) 
#######################################


##SNP in unmapped contigs
#snp.filtered <- snp.filtered[snp.filtered$chrom < 8,] #Keep SNPs that are mapped to the 7 chromosomes
#No need to do this

##Assumes data for a single row (and assumes that all cases are biallelic, missing data is "N")



### * Constructing genetic maps with ASMap


library(ASMap)

testd <- mstmap(mapDHf, dist.fun = "kosambi", trace = TRUE, as.cross = TRUE)

nmar(testd)

chrlen(testd)
       
profileGen(testd, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id =
"Genotype", xo.lambda = 25, layout = c(1, 3), lty = 2)

### ** Testing the tutorial
data(mapBCu) #Load example data

plotMissing(mapBCu) #Plot missing genotypes

sg <- statGen(mapBCu, bychr = FALSE, stat.type = "miss")
mapBC1 <- subset(mapBCu, ind = sg$miss < 1600) #Omit samples that have more than 1600 markers missing

#Identify any potential clones (clones affect segregation distortion)
gc <- genClones(mapBC1, tol = 0.95) #Are there groups that share more than 95% of alleles?
gc$cgd

#Remove clones
cgd <- gc$cgd[-c(1, 4, 5), ]
mapBC2 <- fixClones(mapBC1, cgd, consensus = TRUE)
levels(mapBC2$pheno[[1]])[grep("_", levels(mapBC2$pheno[[1]]))]

#Check segregation distortion
profileMark(mapBC2, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

#Dropping distorted markers
mm <- statMark(mapBC2, stat.type = "marker")$marker$AB #Marker proportions
mapBC3 <- drop.markers(mapBC2, c(markernames(mapBC2)[mm > 0.98], markernames(mapBC2)[mm < 0.2])) #drop the distorted markers

#Pull problematic markers and set them aside
mapBC3 <- pullCross(mapBC3, type = "missing", pars = list(miss.thresh = 0.1))
mapBC3 <- pullCross(mapBC3, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
mapBC3 <- pullCross(mapBC3, type = "co.located")
names(mapBC3)

sum(ncol(mapBC3$missing$data), ncol(mapBC3$seg.dist$data),ncol(mapBC3$co.located$data))

#Do the mapping itself
mapBC4 <- mstmap(mapBC3, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12)
chrlen(mapBC4)

#Plot heatmap of LOD scores and RF
heatMap(mapBC4, lmax = 70) #Need to adjust lmax parameter to have uniform heat

##Looking at genetic distances, there seems to be some individuals that have excess number of crossovers
pg <- profileGen(mapBC4, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)

mapBC5 <- subsetCross(mapBC4, ind = !pg$xo.lambda)
mapBC6 <- mstmap(mapBC5, bychr = TRUE, dist.fun = "kosambi", trace = TRUE, p.value = 1e-12)
chrlen(mapBC6)

#Reasonable number of double cross overs between markers
profileMark(mapBC6, stat.type = c("seg.dist", "prop", "dxo", "recomb"), layout = c(1, 5), type = "l")

#Pushing missing markers back 
mapBC6 <- pushCross(mapBC6, type = "missing", pars = list(miss.thresh = 0.22, max.rf = 0.3))

heatMap(mapBC6, chr = c("L.3", "L.5", "L.8", "L.9"), lmax = 70) #There are genuine linkages betweem groups 3 and 5 and 8 and 9

#Merging linkage groups
mapBC6 <- mergeCross(mapBC6, merge = list(L.3 = c("L.3", "L.5"), L.8 = c("L.8", "L.9")))
names(mapBC6$geno) <- paste("L.", 1:7, sep = "")
mapBC7 <- mstmap(mapBC6, bychr = TRUE, trace = TRUE, dist.fun = "kosambi", p.value = 2)
chrlen(mapBC7)

heatMap(mapBC7, lmax = 70)

##Checking again
pg1 <- profileGen(mapBC7, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)
#Some problems individuals

mapBC8 <- subsetCross(mapBC7, ind = !pg1$xo.lambda)
mapBC9 <- mstmap(mapBC8, bychr = TRUE, dist.fun = "kosambi", trace = TRUE, p.value = 2)
chrlen(mapBC9)

##Check again
profileMark(mapBC9, stat.type = c("seg.dist", "prop", "dxo", "recomb"), layout = c(1, 5), type = "l")

##Pushing back segregation distorted markers
dm <- markernames(mapBC9, "L.2")[statMark(mapBC9, chr = "L.2", stat.type = "marker")$marker$neglog10P > 6] #Remove one marker which looks like non-biological distortion
mapBC10 <- drop.markers(mapBC9, dm)
mapBC11 <- pushCross(mapBC10, type = "seg.distortion", pars = list(seg.ratio = "70:30"))
mapBC12 <- mstmap(mapBC11, bychr = TRUE, trace = TRUE, dist.fun = "kosambi", p.value = 2)
round(chrlen(mapBC12) - chrlen(mapBC9), 5)

mapBC <- pushCross(mapBC12, type = "co.located")
names(mapBC)

plot(mapBC)

plotGeno(mapBC, chr = "L.7", ind = 1:90) #Plots individual genotypes, showing recombination place etc.

### ** Genetic map for Family A
test <- mstmap(famA.asmap, dist.dun = "kosambi", trace = TRUE, as.cross = TRUE)
plotMissing(test)

sg <- statGen(test, bychr = FALSE, stat.type = "miss")
famA.1 <- subset(test, ind = sg$miss < 6000) #Omit samples that have more than 6000 markers missing

#Identify any potential clones (clones affect segregation distortion)
gc <- genClones(famA.1, tol = 0.95) #Are there groups that share more than 95% of alleles?
gc$cgd

#No clones detected

famA.1 <- test

#Check segregation distortion
profileMark(famA.1, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

#Pull problematic markers and set them aside
famA.1 <- pullCross(famA.1, type = "missing", pars = list(miss.thresh = 0.1))
famA.1 <- pullCross(famA.1, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
famA.1 <- pullCross(famA.1, type = "co.located")
names(famA.1)

sum(ncol(famA.1$missing$data), ncol(famA.1$seg.dist$data),ncol(famA.1$co.located$data))

#Do the mapping itself
famA.2 <- mstmap(famA.1, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-6)
chrlen(famA.2)

#Plot heatmap of LOD scores and RF
heatMap(famA.2, lmax = 20) #Need to adjust lmax parameter to have uniform heat (20 seems OK?)

##There seems to be quite reasonable number of double CO's
pg <- profileGen(famA.2, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)

#Reasonable number of double cross overs between markers
pdf(file = "famAdxo.pdf", height = 7, width = 7*1.618)
profileMark(famA.2, stat.type = c("seg.dist", "prop", "dxo", "recomb"), layout = c(1, 5), type = "l")
dev.off()
    
famA.3 <- mstmap(famA.1, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-6)
chrlen(famA.3)

#Plot heatmap of LOD scores and RF
heatMap(famA.3, lmax = 20)


#Checking inds 
plotGeno(famA.2, chr = "L.5", ind = 1:10)

#Viewing the genotypes
view.pgenot(snpdataUG.filt, genotypesUG.filt, ind = "A13", chr = 1)


### *** Plot physical vs. genetic map
famA.3$geno$map

head(famA.3$geno$L.1$data)


##Markers for linkage group 1
markers <- rev(colnames(famA.3$geno$L.1$data) )
m.index <- as.numeric(gsub("marker", "", markers))
physical.pos <- snpdataUG.filt[m.index,1:5]$POS #Chromosome 4
marker.pos <- famA.3$geno$L.1$map #Markers needs to be reversed

##Markers for linkage group 2
markers <- rev(colnames(famA.3$geno$L.2$data))
m.index <- as.numeric(gsub("marker", "", markers))
physical.pos <- snpdataUG.filt[m.index,1:5]$POS #Chromosome 6
marker.pos <- famA.3$geno$L.2$map #Markers needs to be reversed

##Markers for linkage group 3
markers <- colnames(famA.3$geno$L.3$data)
m.index <- as.numeric(gsub("marker", "", markers)) 
physical.pos <- snpdataUG.filt[m.index,1:5]$POS #Chromosome 7
marker.pos <- famA.3$geno$L.3$map #Markers needs to be reversed

##Markers for linkage group 4
markers <- colnames(famA.3$geno$L.4$data)
m.index <- as.numeric(gsub("marker", "", markers)) #Chromosome 5
physical.pos <- snpdataUG.filt[m.index,1:5]$POS
marker.pos <- famA.3$geno$L.4$map #Marker order looks quite OK

##Markers for linkage group 5
markers <- colnames(famA.3$geno$L.5$data)
m.index <- as.numeric(gsub("marker", "", markers)) #Chromosome 1
physical.pos <- snpdataUG.filt[m.index,1:5]$POS
marker.pos <- famA.3$geno$L.5$map #Marker order looks quite OK

##Markers for linkage group 6
markers <- colnames(famA.3$geno$L.6$data)
m.index <- as.numeric(gsub("marker", "", markers)) #Chromosome 3
physical.pos <- snpdataUG.filt[m.index,1:5]$POS
marker.pos <- famA.3$geno$L.6$map #Marker order needs to be reversed

##Markers for linkage group 7
markers <- colnames(famA.3$geno$L.7$data)
m.index <- as.numeric(gsub("marker", "", markers)) #Chromosome 2
physical.pos <- snpdataUG.filt[m.index,1:5]$POS
marker.pos <- famA.3$geno$L.7$map #Marker order needs to be reversed

plot(marker.pos, physical.pos, type ="p", xlab = "Map position (cM)", ylab = "Physical position")

#pdf(file = "FamAmaps.pdf", width = 8*1.618, height = 8) #1.618 is the golden ratio
plot.marker.positions(famA.3, snpdataUG.filt)
#dev.off()

##Need to adjust linkage groups and chromosomes manually
plot.marker.positions <- function(markerdata, snpdata) {
    markers4 <- rev(colnames(markerdata$geno$L.1$data)) #Chromosome 4
    m.index4 <- as.numeric(gsub("marker", "", markers4))
    physical.pos4 <- snpdata[m.index4,]$POS #Chromosome 4
    marker.pos4 <- markerdata$geno$L.1$map #Markers needs to be reversed

    ##Markers for linkage group 2
    markers6 <- rev(colnames(markerdata$geno$L.2$data))
    m.index6 <- as.numeric(gsub("marker", "", markers6))
    physical.pos6 <- snpdata[m.index6,]$POS #Chromosome 6
    marker.pos6 <- markerdata$geno$L.2$map #Markers needs to be reversed

    ##Markers for linkage group 3
    markers7 <- rev(colnames(markerdata$geno$L.3$data))
    m.index7 <- as.numeric(gsub("marker", "", markers7)) 
    physical.pos7 <- snpdata[m.index7,]$POS #Chromosome 7
    marker.pos7 <- markerdata$geno$L.3$map #Markers needs to be reversed

    ##Markers for linkage group 4
    markers5 <- colnames(markerdata$geno$L.4$data)
    m.index5 <- as.numeric(gsub("marker", "", markers5)) #Chromosome 5
    physical.pos5 <- snpdata[m.index5,]$POS
    marker.pos5 <- markerdata$geno$L.4$map #Marker order looks quite OK

    ##Markers for linkage group 5
    markers1 <- colnames(markerdata$geno$L.5$data)
    m.index1 <- as.numeric(gsub("marker", "", markers1)) #Chromosome 1
    physical.pos1 <- snpdata[m.index1,]$POS
    marker.pos1 <- markerdata$geno$L.5$map #Marker order looks quite OK

    ##Markers for linkage group 6
    markers3 <- rev(colnames(markerdata$geno$L.6$data))
    m.index3 <- as.numeric(gsub("marker", "", markers3)) #Chromosome 3
    physical.pos3 <- snpdata[m.index3,]$POS
    marker.pos3 <- markerdata$geno$L.6$map #Marker order needs to be reversed

    ##Markers for linkage group 7
    markers2 <- rev(colnames(markerdata$geno$L.7$data))
    m.index2 <- as.numeric(gsub("marker", "", markers2)) #Chromosome 2
    physical.pos2 <- snpdata[m.index2,]$POS
    marker.pos2 <- markerdata$geno$L.7$map #Marker order needs to be reversed

    ##Combine physical and marker positions
    physical.pos2 <- physical.pos2 + tail(physical.pos1, n = 1)
    physical.pos3 <- physical.pos3 + tail(physical.pos2, n = 1)
    physical.pos4 <- physical.pos4 + tail(physical.pos3, n = 1)
    physical.pos5 <- physical.pos5 + tail(physical.pos4, n = 1)
    physical.pos6 <- physical.pos6 + tail(physical.pos5, n = 1)
    physical.pos7 <- physical.pos7 + tail(physical.pos6, n = 1)

    physical.pos <- c(physical.pos1, physical.pos2, physical.pos3, physical.pos4, physical.pos5, physical.pos6, physical.pos7)

    marker.pos2 <- marker.pos2 + tail(marker.pos1, n = 1)
    marker.pos3 <- marker.pos3 + tail(marker.pos2, n = 1)
    marker.pos4 <- marker.pos4 + tail(marker.pos3, n = 1)
    marker.pos5 <- marker.pos5 + tail(marker.pos4, n = 1)
    marker.pos6 <- marker.pos6 + tail(marker.pos5, n = 1)
    marker.pos7 <- marker.pos7 + tail(marker.pos6, n = 1)

    marker.pos <- c(marker.pos1, marker.pos2, marker.pos3, marker.pos4, marker.pos5, marker.pos6, marker.pos7)

    markerdata <- data.frame(marker.pos, physical.pos)

    ppos <- data.frame(x = c(0, marker.pos2[1], marker.pos2[1], tail(marker.pos2,n=1), tail(marker.pos2,n=1),0, 0, marker.pos4[1], marker.pos4[1], tail(marker.pos4,n=1), tail(marker.pos4,n=1),0, 0, marker.pos6[1], marker.pos6[1], tail(marker.pos6,n=1), tail(marker.pos6,n=1), 0), y = c(physical.pos2[1], physical.pos2[1], 0, 0, tail(physical.pos2,n=1), tail(physical.pos2,n=1), physical.pos4[1], physical.pos4[1], 0, 0, tail(physical.pos4,n=1), tail(physical.pos4,n=1), physical.pos6[1], physical.pos6[1], 0, 0, tail(physical.pos6,n=1), tail(physical.pos6,n=1)), id = rep(factor(c("chr2", "chr4", "chr6")), each = 6))
    
    ggplot(markerdata, aes(x = marker.pos, y = physical.pos)) +
        geom_polygon(aes(x = x, y = y), data = ppos, alpha = 0.2) +
        geom_point(shape = 3) +
        xlab("Map position (cM)") +
        ylab("Physical position") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        annotate("text", x = c(50,50,50,50,50,50,50), y = c(physical.pos2[1]/2, (tail(physical.pos2,n=1) - physical.pos2[1])/2 + physical.pos2[1], (tail(physical.pos3,n=1) - physical.pos3[1])/2 + physical.pos3[1], (tail(physical.pos4,n=1) - physical.pos4[1])/2 + physical.pos4[1], (tail(physical.pos5,n=1) - physical.pos5[1])/2 + physical.pos5[1], (tail(physical.pos6,n=1) - physical.pos6[1])/2 + physical.pos6[1], (tail(physical.pos7,n=1) - physical.pos7[1])/2 + physical.pos7[1]), label = c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr 7"))     
    
    #plot(marker.pos, physical.pos, type ="b", xlab = "Map position (cM)", ylab = "Physical position")

}

    


### ** Genetic map for Family B

test <- mstmap(famB.asmap, dist.dun = "kosambi", trace = TRUE, as.cross = TRUE)
plotMissing(test)

sg <- statGen(test, bychr = FALSE, stat.type = "miss")
famB.1 <- subset(test, ind = sg$miss < 10000) #Omit samples that have more than 1+6000 markers missing

#Identify any potential clones (clones affect segregation distortion)
gc <- genClones(famB.1, tol = 0.95) #Are there groups that share more than 95% of alleles?
gc$cgd
#No clones detected

#Check segregation distortion
profileMark(famB.1, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

#Pull problematic markers and set them aside
famB.1 <- pullCross(famB.1, type = "missing", pars = list(miss.thresh = 0.1))
famB.1 <- pullCross(famB.1, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
famB.1 <- pullCross(famB.1, type = "co.located")
names(famB.1)

sum(ncol(famB.1$missing$data), ncol(famB.1$seg.dist$data),ncol(famB.1$co.located$data))

#Do the mapping itself
famB.2 <- mstmap(famB.1, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-5)
chrlen(famB.2)

#Plot heatmap of LOD scores and RF
heatMap(famB.2, lmax = 10) #Need to adjust lmax parameter to have uniform heat (20 seems OK?)

##There seems to be quite reasonable number of double CO's
pg <- profileGen(famB.2, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)

#Reasonable number of double cross overs between markers
pdf(file = "famBdxo.pdf", height = 7, width = 7*1.618)
profileMark(famB.2, stat.type = c("seg.dist", "prop", "dxo", "recomb"), layout = c(1, 5), type = "l")
dev.off()

#Checking inds 
plotGeno(famB.2, chr = "L.5", ind = 1:10)

#Viewing the genotypes
view.pgenot(snpdataUG.filt, genotypesUG.filt, ind = "B13", chr = 1)

### *** Plot physical vs. genetic position

plot(marker.pos, physical.pos, type ="p", xlab = "Map position (cM)", ylab = "Physical position")


##Need to adjust linkage groups and chromosomes manually
plot.marker.positions <- function(markerdata, snpdata) {
##Markers for linkage group 1
markers3 <- rev(colnames(famB.2$geno$L.1$data) )
m.index3 <- as.numeric(gsub("marker", "", markers3))
physical.pos3 <- snpdataUG.filt[m.index3,1:5]$POS #Chromosome 3
marker.pos3 <- famB.2$geno$L.1$map #Markers needs to be reversed
#
##Markers for linkage group 2
markers4 <- colnames(famB.2$geno$L.2$data)
m.index4 <- as.numeric(gsub("marker", "", markers4))
physical.pos4 <- snpdataUG.filt[m.index4,1:5]$POS #Chromosome 4
marker.pos4 <- famB.2$geno$L.2$map
#
##Markers for linkage group 3
markers1 <- colnames(famB.2$geno$L.3$data)
m.index1 <- as.numeric(gsub("marker", "", markers1)) 
physical.pos1 <- snpdataUG.filt[m.index1,1:5]$POS #Chromosome 1
marker.pos1 <- famB.2$geno$L.3$map
#
##Markers for linkage group 4
markers2 <- rev(colnames(famB.2$geno$L.4$data))
m.index2 <- as.numeric(gsub("marker", "", markers2)) #Chromosome 2
physical.pos2 <- snpdataUG.filt[m.index2,1:5]$POS
marker.pos2 <- famB.2$geno$L.4$map #Marker order needs to be reversed
#
##Markers for linkage group 5
markers5 <- rev(colnames(famB.2$geno$L.5$data))
m.index5 <- as.numeric(gsub("marker", "", markers5)) #Chromosome 5
physical.pos5 <- snpdataUG.filt[m.index5,1:5]$POS
marker.pos5 <- famB.2$geno$L.5$map #Marker order looks quite OK
#
##Markers for linkage group 6
markers6 <- colnames(famB.2$geno$L.6$data)
m.index6 <- as.numeric(gsub("marker", "", markers6)) #Chromosome 6
physical.pos6 <- snpdataUG.filt[m.index6,1:5]$POS
marker.pos6 <- famB.2$geno$L.6$map #Marker order needs to be reversed
#
##Markers for linkage group 7
markers7 <- colnames(famB.2$geno$L.7$data)
m.index7 <- as.numeric(gsub("marker", "", markers7)) #Chromosome 7
physical.pos7 <- snpdataUG.filt[m.index7,1:5]$POS
marker.pos7 <- famB.2$geno$L.7$map #Marker order needs to be reversed
#
    ##Combine physical and marker positions
    physical.pos2 <- physical.pos2 + tail(physical.pos1, n = 1)
    physical.pos3 <- physical.pos3 + tail(physical.pos2, n = 1)
    physical.pos4 <- physical.pos4 + tail(physical.pos3, n = 1)
    physical.pos5 <- physical.pos5 + tail(physical.pos4, n = 1)
    physical.pos6 <- physical.pos6 + tail(physical.pos5, n = 1)
    physical.pos7 <- physical.pos7 + tail(physical.pos6, n = 1)
#
    physical.pos <- c(physical.pos1, physical.pos2, physical.pos3, physical.pos4, physical.pos5, physical.pos6, physical.pos7)
#
    marker.pos2 <- marker.pos2 + tail(marker.pos1, n = 1)
    marker.pos3 <- marker.pos3 + tail(marker.pos2, n = 1)
    marker.pos4 <- marker.pos4 + tail(marker.pos3, n = 1)
    marker.pos5 <- marker.pos5 + tail(marker.pos4, n = 1)
    marker.pos6 <- marker.pos6 + tail(marker.pos5, n = 1)
    marker.pos7 <- marker.pos7 + tail(marker.pos6, n = 1)
#
    marker.pos <- c(marker.pos1, marker.pos2, marker.pos3, marker.pos4, marker.pos5, marker.pos6, marker.pos7)
#
    markerdata <- data.frame(marker.pos, physical.pos)
#
    ppos <- data.frame(x = c(0, marker.pos2[1], marker.pos2[1], tail(marker.pos2,n=1), tail(marker.pos2,n=1),0, 0, marker.pos4[1], marker.pos4[1], tail(marker.pos4,n=1), tail(marker.pos4,n=1),0, 0, marker.pos6[1], marker.pos6[1], tail(marker.pos6,n=1), tail(marker.pos6,n=1), 0), y = c(physical.pos2[1], physical.pos2[1], 0, 0, tail(physical.pos2,n=1), tail(physical.pos2,n=1), physical.pos4[1], physical.pos4[1], 0, 0, tail(physical.pos4,n=1), tail(physical.pos4,n=1), physical.pos6[1], physical.pos6[1], 0, 0, tail(physical.pos6,n=1), tail(physical.pos6,n=1)), id = rep(factor(c("chr2", "chr4", "chr6")), each = 6))
#
    ggplot(markerdata, aes(x = marker.pos, y = physical.pos)) +
        geom_polygon(aes(x = x, y = y), data = ppos, alpha = 0.2) +
        geom_point(shape = 3) +
        xlab("Map position (cM)") +
        ylab("Physical position") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        annotate("text", x = c(25,25,25,25,25,25,25), y = c(physical.pos2[1]/2, (tail(physical.pos2,n=1) - physical.pos2[1])/2 + physical.pos2[1], (tail(physical.pos3,n=1) - physical.pos3[1])/2 + physical.pos3[1], (tail(physical.pos4,n=1) - physical.pos4[1])/2 + physical.pos4[1], (tail(physical.pos5,n=1) - physical.pos5[1])/2 + physical.pos5[1], (tail(physical.pos6,n=1) - physical.pos6[1])/2 + physical.pos6[1], (tail(physical.pos7,n=1) - physical.pos7[1])/2 + physical.pos7[1]), label = c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr 7"))     
    #
    #plot(marker.pos, physical.pos, type ="b", xlab = "Map position (cM)", ylab = "Physical position")
#
}

pdf(file = "FamBmaps.pdf", width = 8*1.618, height = 8) #1.618 is the golden ratio
plot.marker.positions(famB.2, snpdataUG.filt)
dev.off()



### ** Genetic map for Family C

test <- mstmap(famC.asmap, dist.dun = "kosambi", trace = TRUE, as.cross = TRUE)
plotMissing(test)

sg <- statGen(test, bychr = FALSE, stat.type = "miss")
famC.1 <- subset(test, ind = sg$miss < 10000) #Omit samples that have more than 10000 markers missing

#Identify any potential clones (clones affect segregation distortion)
gc <- genClones(famC.1, tol = 0.95) #Are there groups that share more than 95% of alleles?
gc$cgd
#No clones detected

#Check segregation distortion
profileMark(famC.1, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

#Pull problematic markers and set them aside
famC.1 <- pullCross(famC.1, type = "missing", pars = list(miss.thresh = 0.1))
famC.1 <- pullCross(famC.1, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
famC.1 <- pullCross(famC.1, type = "co.located")
names(famC.1)

sum(ncol(famC.1$missing$data), ncol(famC.1$seg.dist$data),ncol(famC.1$co.located$data))

#Do the mapping itself
famC.2 <- mstmap(famC.1, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-6)
chrlen(famC.2)

#Plot heatmap of LOD scores and RF
heatMap(famC.2, lmax = 10) #Need to adjust lmax parameter to have uniform heat (20 seems OK?)

##There seems to be quite reasonable number of double CO's
pg <- profileGen(famC.2, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)

#Reasonable number of double cross overs between markers
pdf(file = "famCdxo.pdf", height = 7, width = 7*1.618)
profileMark(famC.2, stat.type = c("seg.dist", "prop", "dxo", "recomb"), layout = c(1, 5), type = "l")
dev.off()

#Checking inds 
plotGeno(famC.2, chr = "L.5", ind = 1:10)

#Viewing the genotypes
view.pgenot(snpdataUG.filt, genotypesUG.filt, ind = "C20", chr = 7)


### ** Genetic map for Family D

test <- mstmap(famD.asmap, dist.dun = "kosambi", trace = TRUE, as.cross = TRUE)
plotMissing(test)

sg <- statGen(test, bychr = FALSE, stat.type = "miss")
famD.1 <- subset(test, ind = sg$miss < 10000) #Omit samples that have more than 1+6000 markers missing

#Identify any potential clones (clones affect segregation distortion)
gc <- genClones(famD.1, tol = 0.95) #Are there groups that share more than 95% of alleles?
gc$cgd
#No clones detected

#Check segregation distortion
profileMark(famD.1, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

#Pull problematic markers and set them aside
famD.1 <- pullCross(famD.1, type = "missing", pars = list(miss.thresh = 0.1))
famD.1 <- pullCross(famD.1, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
famD.1 <- pullCross(famD.1, type = "co.located")
names(famD.1)

sum(ncol(famD.1$missing$data), ncol(famD.1$seg.dist$data),ncol(famD.1$co.located$data))

#Do the mapping itself
famD.2 <- mstmap(famD.1, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-6)
chrlen(famD.2)

#Plot heatmap of LOD scores and RF
heatMap(famD.2, lmax = 10) #Need to adjust lmax parameter to have uniform heat (20 seems OK?)

##There seems to be quite reasonable number of double CO's
pg <- profileGen(famD.2, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)

#Reasonable number of double cross overs between markers
pdf(file = "famDdxo.pdf", height = 7, width = 7*1.618)
profileMark(famD.2, stat.type = c("seg.dist", "prop", "dxo", "recomb"), layout = c(1, 5), type = "l")
dev.off()

#Checking inds 
plotGeno(famD.2, chr = "L.5", ind = 1:10)

#Viewing the genotypes
view.pgenot(snpdataUG.filt, genotypesUG.filt, ind = "D13", chr = 1)

### ** Genetic map for Family E

test <- mstmap(famE.asmap, dist.dun = "kosambi", trace = TRUE, as.cross = TRUE)
plotMissing(test)

sg <- statGen(test, bychr = FALSE, stat.type = "miss")
famE.1 <- subset(test, ind = sg$miss < 10000) #Omit samples that have more than 1+6000 markers missing

#Identify any potential clones (clones affect segregation distortion)
gc <- genClones(famE.1, tol = 0.95) #Are there groups that share more than 95% of alleles?
gc$cgd
#No clones detected

#Check segregation distortion
profileMark(famE.1, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

#Pull problematic markers and set them aside
famE.1 <- pullCross(famE.1, type = "missing", pars = list(miss.thresh = 0.1))
famE.1 <- pullCross(famE.1, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
famE.1 <- pullCross(famE.1, type = "co.located")
names(famE.1)

sum(ncol(famE.1$missing$data), ncol(famE.1$seg.dist$data),ncol(famE.1$co.located$data))

#Do the mapping itself
famE.2 <- mstmap(famE.1, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-6)
chrlen(famE.2)

#Plot heatmap of LOD scores and RF
heatMap(famE.2, lmax = 10) #Need to adjust lmax parameter to have uniform heat (20 seems OK?)

##There seems to be quite reasonable number of double CO's
pg <- profileGen(famE.2, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)

#Reasonable number of double cross overs between markers
pdf(file = "famEdxo.pdf", height = 7, width = 7*1.618)
profileMark(famE.2, stat.type = c("seg.dist", "prop", "dxo", "recomb"), layout = c(1, 5), type = "l")
dev.off()

#Checking inds 
plotGeno(famE.2, chr = "L.5", ind = 1:10)

#Viewing the genotypes
view.pgenot(snpdataUG.filt, genotypesUG.filt, ind = "E13", chr = 1)


### ** Genetic map for Family G

test <- mstmap(famG.asmap, dist.dun = "kosambi", trace = TRUE, as.cross = TRUE)
plotMissing(test)

sg <- statGen(test, bychr = FALSE, stat.type = "miss")
famG.1 <- subset(test, ind = sg$miss < 10000) #Omit samples that have more than 10000 markers missing

#Identify any potential clones (clones affect segregation distortion)
gc <- genClones(famG.1, tol = 0.95) #Are there groups that share more than 95% of alleles?
gc$cgd
#Samples G35 and G29 are potential clones
#Samples G64 and G58 are potential clones
subind <- sg$miss < 10000
subind[c(21,52)] <- FALSE
famG.1 <- subset(test, ind = subind)

#Check segregation distortion
profileMark(famG.1, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

#Pull problematic markers and set them aside
famG.1 <- pullCross(famG.1, type = "missing", pars = list(miss.thresh = 0.1))
famG.1 <- pullCross(famG.1, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
famG.1 <- pullCross(famG.1, type = "co.located")
names(famG.1)

sum(ncol(famG.1$missing$data), ncol(famG.1$seg.dist$data),ncol(famG.1$co.located$data))

#Do the mapping itself
famG.2 <- mstmap(famG.1, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-6)
chrlen(famG.2)

#Plot heatmap of LOD scores and RF
heatMap(famG.2, lmax = 10) #Need to adjust lmax parameter to have uniform heat (20 seems OK?)

##There seems to be quite an excess amount of double CO's
pg <- profileGen(famG.2, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)

#Reasonable number of double cross overs between markers
pdf(file = "famGdxo.pdf", height = 7, width = 7*1.618)
profileMark(famG.2, stat.type = c("seg.dist", "prop", "dxo", "recomb"), layout = c(1, 5), type = "l")
dev.off()

#Checking inds 
plotGeno(famG.2, chr = "L.5", ind = 1:10)

#Viewing the genotypes
view.pgenot(snpdataUG.filt, genotypesUG.filt, ind = "G62", chr = 7)

#Need to do some manual curation
pg$stat$dxo #This shows the number of couble xo's for each individual


