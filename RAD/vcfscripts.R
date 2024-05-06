##These scripts deal with vcf files and genotypes

#Used in dealing with RAD-seq data

##Process .vcf genotypes
##Input is the .vcf data, and cols are the genotype columns, i.e. 11:100
process.vcf.genotypes <- function(data,cols) {
    genotypes <- data[,cols]
    genotypes <- apply(genotypes, 2, list) #This makes a list of lists
    genotypes <- lapply(genotypes, unlist) #This recudes it to a single list
    genotypes <- lapply(genotypes, matrix, ncol = 1)
    for(i in 1:length(genotypes)) {colnames(genotypes[[i]]) <- "dat"} #change colnames
    genotypes <- lapply(genotypes, data.frame)
    genotypes <- lapply(genotypes, separate, col = "dat", into = c("GT", "DP", "AD", "QL", "GL"), sep = ":")
    return(genotypes)
}

##Version for unified genotyper
process.vcf.genotypesUG <- function(data,cols) {
    genotypes <- data[,cols]
    genotypes <- apply(genotypes, 2, list) #This makes a list of lists
    genotypes <- lapply(genotypes, unlist) #This recudes it to a single list
    genotypes <- lapply(genotypes, matrix, ncol = 1)
    for(i in 1:length(genotypes)) {colnames(genotypes[[i]]) <- "dat"} #change colnames
    genotypes <- lapply(genotypes, data.frame)
    genotypes <- lapply(genotypes, separate, col = "dat", into = c("GT", "AD", "DP", "GQ", "PL"), sep = ":")
    return(genotypes)
}

##Calculate average genotype quality and sequencing depth for each lopcus across all individuals
avg.quality.depth.vcf <- function(genotypes) {
    nloc <- dim(genotypes[[1]])[1] #count number of loci
    res.mat <- matrix(rep(0,nloc*2), ncol = 2)
    colnames(res.mat) <- c("avgGQ", "avgDP")

    GQs <- matrix(rep(0, nloc*length(genotypes)), ncol = length(genotypes))
    DPs <- matrix(rep(0, nloc*length(genotypes)), ncol = length(genotypes))
    
    for(i in 1:nloc) { #
        for(j in 1:length(genotypes)) {
                GQs[i,j] <- as.numeric(genotypes[[j]][i,4]) #store genotype quality
                DPs[i,j] <- as.numeric(genotypes[[j]][i,2]) #store allele depth
            }
    }

    #Count means
    res.mat[,1] <- apply(GQs, 1, mean, na.rm = TRUE)
    res.mat[,2] <- apply(DPs, 1, mean, na.rm = TRUE)

    return(res.mat)
}

##Calculate average genotype quality and sequencing depth for each lopcus across all individuals (Unified genotyper
avg.quality.depth.vcf.UG <- function(genotypes) {
    nloc <- dim(genotypes[[1]])[1] #count number of loci
    res.mat <- matrix(rep(0,nloc*2), ncol = 2)
    colnames(res.mat) <- c("avgGQ", "avgDP")

    GQs <- matrix(rep(0, nloc*length(genotypes)), ncol = length(genotypes))
    DPs <- matrix(rep(0, nloc*length(genotypes)), ncol = length(genotypes))
    
    for(i in 1:nloc) { #
        for(j in 1:length(genotypes)) {
                GQs[i,j] <- as.numeric(genotypes[[j]][i,4]) #store genotype quality
                DPs[i,j] <- as.numeric(genotypes[[j]][i,3]) #store allele depth
            }
    }

    #Count means
    res.mat[,1] <- apply(GQs, 1, mean, na.rm = TRUE)
    res.mat[,2] <- apply(DPs, 1, mean, na.rm = TRUE)

    return(res.mat)
}


##Check the proportion of hets
hets.vcf <- function(genotypes) {
    nloc <- dim(genotypes[[1]])[1] #count number of loci
    res.mat <- matrix(rep(0,nloc*1), ncol = 1)
    colnames(res.mat) <- c("propHET")

    GTs <- matrix(rep(0, nloc*length(genotypes)), ncol = length(genotypes))

    for(i in 1:nloc) { #
        for(j in 1:length(genotypes)) {
                GTs[i,j] <- as.character(genotypes[[j]][i,1]) #store genotype
                
            }
    }

    #Calculate the proportion of heterozygotes of those genotypes that are not missing
    hetprop <- function(x) { sum(x == "0/1") / (length(x) - sum(x == "./.")) }

    res.mat[,1] <- apply(GTs, 1, hetprop)

    return(res.mat)
}

##Check the proportion of hets
hets.vcf.UG <- function(genotypes) {
    nloc <- dim(genotypes[[1]])[1] #count number of loci
    res.mat <- matrix(rep(0,nloc*1), ncol = 1)
    colnames(res.mat) <- c("propHET")

    GTs <- matrix(rep(0, nloc*length(genotypes)), ncol = length(genotypes))

    for(i in 1:nloc) { #
        for(j in 1:length(genotypes)) {
                GTs[i,j] <- as.character(genotypes[[j]][i,1]) #store genotype
                
            }
    }

    #Calculate the proportion of heterozygotes of those genotypes that are not missing
    hetprop <- function(x) {
        y <- strsplit(x, "/")
        x1 <- sapply(y, '[[',1)
        x2 <- sapply(y, '[[',2)
        hets <- x1 != x2
        sum(hets) / (length(x) - sum(x == "./.")) }

    res.mat[,1] <- apply(GTs, 1, hetprop)

    return(res.mat)
}

##This function checks whether a genotype is heterozygous
##Returns TRUE if it is het
is.het.UG <- function(x) {
    y <- strsplit(x, "/")
    x1 <- sapply(y, '[[',1)
    x2 <- sapply(y, '[[',2)
    check <- x1 != x2
    return(check)
}

calc.MAF.UG <- function(x) {
    y <- 1 - x
    res <- cbind(x,y)
    results <- apply(res, MARGIN = 1, min)
    return(results)
}


##This function calculates the average sequencing depths for both genotypes and also calculates the ratio of allelic sequencing depths
allelic.depth.balance <- function(genotypes) {
    nloc <- nrow(genotypes[[1]])
    genotypes <- lapply(genotypes, separate, col = "AD", into = c("AD0", "AD1"), sep = ",")

    #For each locus, pick all allelic depths for each 0/0 genotype, and allelic depth for each 1/1 genotype
    res.mat <- matrix(rep(0, nloc*3), ncol = 3)
    colnames(res.mat) <- c("AD:0/0", "AD:1/1", "ratio")

    #Pick genotypes and allelic depths for each locus
    GTs <- unname(sapply(genotypes, '[[',1))
    AD0s <- unname(sapply(genotypes, '[[',3))
    AD1s <- unname(sapply(genotypes, '[[',4))

    #Calculate mean allelic depths for both genotypes 
    for(i in 1:nloc) {
        res.mat[i,1] <- mean(as.numeric(AD0s[i, GTs[i,] == "0/0"]), na.rm = TRUE)
        res.mat[i,2] <- mean(as.numeric(AD1s[i, GTs[i,] == "1/1"]), na.rm = TRUE)
    }

    res.mat[,3] <- res.mat[,1] / res.mat[,2]

    return(res.mat)
}

##This function calculates the average sequencing depths for both genotypes and also calculates the ratio of allelic sequencing depths
allelic.depth.balance.UG <- function(genotypes) {
    nloc <- nrow(genotypes[[1]])
    genotypes <- lapply(genotypes, separate, col = "AD", into = c("AD0", "AD1", "AD2"), sep = ",")

    #For each locus, pick all allelic depths for each 0/0 genotype, and allelic depth for each 1/1 genotype, and 2/2 genotype
    res.mat <- matrix(rep(0, nloc*4), ncol = 4)
    colnames(res.mat) <- c("AD:0/0", "AD:1/1", "AD:2/2", "ratio")

    #Pick genotypes and allelic depths for each locus
    GTs <- unname(sapply(genotypes, '[[',1))
    AD0s <- unname(sapply(genotypes, '[[',2))
    AD1s <- unname(sapply(genotypes, '[[',3))
    AD2s <- unname(sapply(genotypes, '[[',4))

    #Calculate mean allelic depths for all genotypes 
    for(i in 1:nloc) {
        res.mat[i,1] <- mean(as.numeric(AD0s[i, GTs[i,] == "0/0"]), na.rm = TRUE)
        res.mat[i,2] <- mean(as.numeric(AD1s[i, GTs[i,] == "1/1"]), na.rm = TRUE)
        res.mat[i,3] <- mean(as.numeric(AD2s[i, GTs[i,] == "2/2"]), na.rm = TRUE)

        #Which of the three genotypes are the most common?
        if(sum(is.nan(res.mat[i,1:3])) == 2) { res.mat[i,4] <- 1 }
        if(sum(is.nan(res.mat[i,1:3])) == 1) {
            ind <- which(is.nan(res.mat[i,1:3]))
            twoleft <- res.mat[i, -c(ind,4)]
            res.mat[i,4] <- twoleft[1] / twoleft[2] }
        if(sum(is.nan(res.mat[i,1:3])) == 0) {
            ind <- which(min(res.mat[i,1:3]) == res.mat[i,1:3]) #This is the rarest to drop
            twoleft <- res.mat[i, -c(ind,4)]
            res.mat[i,4] <- twoleft[1] / twoleft[2] }
        
    }

   return(res.mat)
} ###This still has some problem! Need to debug!



##This function converts genotype calls (1/1, 0/0) to SNP calls based on ref and alt bases
##At the moment is deals only with haploid genotype calls (Neurospora)
vcf.gts2snps <- function(snpdata, genotypes) {

    ##Build a data.frame for subsetting
    index <- 1:nrow(snpdata)
    bases <- data.frame(index, snpdata$REF, snpdata$ALT)
    colnames(bases)[2:3] <- c("REF", "ALT")
    
    #loop over all individuals
    for(j in 1:length(genotypes)) {
        refsnp <- genotypes[[j]][,1] == "0/0"
        altsnp <- genotypes[[j]][,1] == "1/1"
        nasnp <- genotypes[[j]][,1] == "0/1" | genotypes[[j]][,1] == "./."
        selectref <- bases[refsnp,-3] #Drop alt column
        selectalt <- bases[altsnp,-2] #Drop ref column
        colnames(selectref)[2] <- "SNP"
        colnames(selectalt)[2] <- "SNP"
        #Check are there any heterozygotes in the first place, combine bases and sort
        if(any(nasnp) == TRUE) {
        selectna <- cbind(bases[nasnp,1], rep(NA, length(bases[nasnp,1]))) #
        colnames(selectna) <- c("index", "SNP")
        SNP <- arrange(rbind(selectref, selectalt, selectna), index)[,2] } else {
            SNP <- arrange(rbind(selectref, selectalt), index)[,2] }

        #Store new data
        genotypes[[j]] <- cbind(genotypes[[j]], SNP)
        genotypes[[j]]$SNP <- as.character(genotypes[[j]]$SNP)
    }

    return(genotypes)
}

##New version of vcf.gts2snps to deal with UG formats
##Assumes that ALT columns has max two alleles (so max two alternative alleles segregating
vcf.gts2snps.UG <- function(snpdata, genotypes) {

    ##Build a data.frame
    index <- 1:nrow(snpdata)
    snpdata <- separate(snpdata, col = ALT, into = c("ALT1", "ALT2"), sep = ",")
    bases <- data.frame(index, snpdata$REF, snpdata$ALT1, snpdata$ALT2) ###
    colnames(bases)[2:4] <- c("REF", "ALT1", "ALT2")

    #Loop over all individuals
    for( j in 1:length(genotypes)) {
        refsnp <- genotypes[[j]][,1] == "0/0"
        altsnp1 <- genotypes[[j]][,1] == "1/1"
        altsnp2 <- genotypes[[j]][,1] == "2/2"
        nasnp <- genotypes[[j]][,1] == "0/1" | genotypes[[j]][,1] == "./." | genotypes[[j]][,1] == "0/2" | genotypes[[j]][,1] == "1/2"
        selectref <- bases[refsnp,-c(3,4)] #Drop alt columns
        selectalt1 <- bases[altsnp1,-c(2,4)] #Drop ref column and alt2
        selectalt2 <- bases[altsnp2,-c(2,3)] #Drop ref column and alt1
        colnames(selectref)[2] <- "SNP"
        colnames(selectalt1)[2] <- "SNP"
        colnames(selectalt2)[2] <- "SNP"
        #Check are there any heterozygotes in the first place, combine bases and sort
        if(any(nasnp) == TRUE) {
        selectna <- cbind(bases[nasnp,1], rep(NA, length(bases[nasnp,1]))) #
        colnames(selectna) <- c("index", "SNP")
        SNP <- arrange(rbind(selectref, selectalt1, selectalt2, selectna), index)[,2] } else {
            SNP <- arrange(rbind(selectref, selectalt1, selectalt2), index)[,2] }

        #Store new data
        genotypes[[j]] <- cbind(genotypes[[j]], SNP)
        genotypes[[j]]$SNP <- as.character(genotypes[[j]]$SNP)
    }
    
    return(genotypes)
}
        
    

##This function converts SNP calls based on parent genotypes into R/qtl input format
##i.e. 1 and 2
vcf.parentgenot <- function(genotypes,parent1ind,parent2ind) {

    print(paste("Parent 'A' is", names(genotypes)[parent1ind], "and parent 'B' is", names(genotypes)[parent2ind]))  
    index <- 1:nrow(genotypes[[parent1ind]])
    #parent1 <- rep("A", nrow(genotypes[[parent1ind]]))
    #parent2 <- rep("B", nrow(genotypes[[parent2ind]]))

    parent1snp <- as.character(genotypes[[parent1ind]]$SNP)
    parent1snp <- cbind(index, parent1snp)
    parent2snp <- as.character(genotypes[[parent2ind]]$SNP)
    parent2snp <- cbind(index, parent2snp)

    #Expects an input vector of size 2
    other.parent.missing.check <- function(input) {
        if(all(is.na(input)) == TRUE) { input[1:2] <- c(FALSE, FALSE) } #If sample is missing
        if(is.na(input[1]) == TRUE & input[2] == TRUE) { input[1] <- FALSE }
        if(is.na(input[2]) == TRUE & input[1] == TRUE) { input[2] <- FALSE }
        if(is.na(input[1]) == TRUE & input[2] == FALSE) { input[1] <- TRUE }
        if(is.na(input[2]) == TRUE & input[1] == FALSE) { input[2] <- TRUE }
        return(input)
    }
    
    
    for(j in 1:length(genotypes)) {
        GENOP <- rep(0, length(index))
        checkp1 <- genotypes[[j]]$SNP == parent1snp[,2]
        checkp2 <- genotypes[[j]]$SNP == parent2snp[,2]
                
        #Check for cases where other parent was missing
        checkcomb <- cbind(checkp1, checkp2)
        testcomb <- t(apply(checkcomb, 1, other.parent.missing.check))
                
        #Fix cases where sample information was missing
        checkna <- is.na(genotypes[[j]]$SNP)
        
        #Make final genotype vector
        GENOP[testcomb[,1]] <- "A" #parent 1 gets genotype A
        GENOP[testcomb[,2]] <- "B" #parent 2 gets genotype B
        GENOP[checkna] <- NA

        #Store final genotype vector
        genotypes[[j]] <- cbind(genotypes[[j]], GENOP)
        genotypes[[j]]$GENOP <- as.character(genotypes[[j]]$GENOP)
    }

    return(genotypes)
}


##This function exports vcf type data into hapmap
export.vcf2hapmap <- function(snpdata, genotypes) {

    nloc <- nrow(snpdata) #Number of loci (numbers of rows)
    nind <- length(genotypes) #Number of samples (needs this + 11 columns)
    res.mat <- data.frame(matrix(rep(0, nloc*(nind+11)), ncol = nind + 11))

    namesforcols <- c("rs", "alleles", "chrom", "pos", "strand", "assembly", "center", "protLSID", "assayLSID", "panel", "QCcode", names(genotypes))
    colnames(res.mat) <- namesforcols #Set column names

    res.mat[,1] <- paste("locus", 1:nloc, sep = "") #SNP names
    res.mat[,c(2,5,7:11)] <- NA #Missing values for many columns
    res.mat[,6] <- "NC12" #Set the genome assembly

    #Setting chromosome numbers and genome positions
    if(length(unique(snpdata[,1])) != 7) { stop("Chromosome number different from 7, check input data!") }
    res.mat[,3] <- gsub("Supercontig_12.", "", snpdata[,1]) #Store chromosome numbers
    res.mat[,4] <- snpdata[,2] #Store genome position

    for(j in 1:nind) { res.mat[,j+11] <- genotypes[[j]]$SNP } #Loop over inds and store genotypes

    return(res.mat)

}

#This function exports parent genotypes to ASMap format
#Assumes that parent genotypes are in "A" and "B" and that last two individuals in the genotype list are the parents
export.vcf2asmap <- function(snpdata, genotypes) {
    
    nloc <- nrow(snpdata) #Number of loci (numbers of rows)
    nind <- length(genotypes) - 2 #Number of samples but parents are dropped
    res.mat <- data.frame(matrix(rep(0, nloc*nind), ncol = nind))

    row.names(res.mat) <- paste("marker", 1:nloc, sep = "") #SNP names
    names(res.mat) <- names(genotypes)[1:nind]

    for(j in 1:nind) {
        res.mat[,j] <- genotypes[[j]]$GENOP
        res.mat[is.na(res.mat[,j]),j] <- "U" #"U" is missing data in ASMap...
    }

    return(res.mat)
}

##Unified Genotyper gives a bit different .vcf output
##Count the number of samples with data per locus
NS.count <- function(aineisto) { sum(aineisto != "./.") }

##This function sets genotypes with low sequencing depth to NA
##Argument min.seq.depth controls the minimum sequencing depth required for a genotype to be called
##default = 6, so genotypes with <= 5 reads are set to NA
quality.control.genot <- function(genotypes, min.seq.depth = 6) {
    ##Loop over all individuals
    for(i in 1:length(genotypes)) {
            index <- as.numeric(genotypes[[i]]$DP) < min.seq.depth #Check for depth smaller than min.seq.depth
            genotypes[[i]]$GT[index] <- "./." #Replace genotypes with missing values
        }

        return(genotypes)
}

#ccounts <- str_count(snpdata$ALT, ",") ##If there is more than one comma there are multiple alleles

##This function checks for multiple segragating alleles
##Returns a logical vector of loci with more than two segregating alleles
check.allele.numbers <- function(genotypes) {
    GTs <- unname(sapply(genotypes, '[[',1))
    segregating.genot <- function(x) { length(table(x, exclude = c("./.", "0/1", "0/2", "1/2"))) }

    n.alleles <- apply(GTs, MARGIN = 1, segregating.genot) #Number of segregating alleles
    is.monomorph <- n.alleles == 1
    is.not.biall <- n.alleles > 2
    results <- data.frame(n.alleles, is.monomorph, is.not.biall)
    return(results)
}

##This function checks if locus is an indel or SNP
is.indel <- function(snpdata) {
    refchar <- nchar(snpdata$REF)
    altchar <- nchar(snpdata$ALT)
    reftest <- refchar > 1 #If ref is longer than 1 -> indel
    #Need to take into account a situation where alt has two SNPs that are both different from ref, i.e ALT = A,C
    alttest <- altchar == 2 | altchar > 3 #If alt = 2 | alt > 3 -> indel
    alttest2 <- altchar == 3 & grepl(",", snpdata$ALT) == FALSE #indel of size 3
    indel <- as.logical(reftest + alttest + alttest2)
    return(indel)
}



#This function is needed in inferring the chromosome segments from GENOP genotypes
#Make a function "find.segments" out of this
find.segments <- function(parentgenos, genotype) {
    isA <- parentgenos == genotype
    seg.res <- c(0,0) #initialize
    segment.l <- FALSE #initialize, logical indicating are we in a segment or not?
    for(i in 1:length(isA)) {
        if(isA[i] == T & segment.l == FALSE & is.na(isA[i]) == FALSE) { #start segment if not in segment
            segment.start <- i
            segment.l <- TRUE
            segment.index <- i
        }
                                        #if(isA[i] == T & i == 1 & is.na(isA[i]) == FALSE) { #This can be simplifed, drop this condition
                                        #    segment.start <- i
                                        #    segment.l <- TRUE
                                        #}
        if(isA[i] == T & i == length(isA) & is.na(isA[i]) == FALSE) { #End segment if end of chr
            segment.end <- i
            seg.res <- rbind(seg.res, c(segment.start, segment.end))
            segment.l <- FALSE
        }
                                        #if(isA[i] == T & i > 1 & segment.l == FALSE & is.na(isA[i]) == FALSE) {
                                        #    segment.start <- i
                                        #    segment.l <- TRUE
                                        #}
        if(isA[i] == F & segment.l == TRUE & is.na(isA[i]) == FALSE) { #End segment if change
            segment.end <- segment.index
            seg.res <- rbind(seg.res, c(segment.start, segment.end))
            segment.l <- FALSE
        }
#
        if(segment.l == TRUE & isA[i] == T & is.na(isA[i]) == FALSE) { #Coordinate of current segment
            segment.index <- i
        }
    }
    if(all(seg.res == 0) == FALSE) {
    seg.res <- seg.res[-1,,drop = FALSE] #Drop initial row (keeps as matrix even with 1 row)
} else { seg.res <- matrix(seg.res, ncol = 2) }
    return(seg.res)
}


##This function infers genomic coordinates of recombination breakpoints
##Note: I assume that segments are done chromosome by chromosome
segments.2.genomic.coordinates <- function(segmentcoordA, segmentcoordB, snpdata, chr) {

    nsegmentA <- nrow(segmentcoordA) #Number of segments
    genomicA <- matrix(rep(0,nsegmentA*2), ncol = 2)

    nsegmentB <- nrow(segmentcoordB) #Number of segments
    genomicB <- matrix(rep(0,nsegmentB*2), ncol = 2)

    #Check that segments exist
    if(nrow(segmentcoordA) == 1 & all(segmentcoordA == 0)) {
        noA <- T } else { noA <- F }

    if(nrow(segmentcoordB) == 1 & all(segmentcoordB == 0)) {
        noB <- T } else { noB <- F }

    #Chromosome ends are treated differently
    ends <- matrix(c(1, 9798894, 2, 4478684, 3, 5274803, 4, 6000762, 5, 6436247, 6, 4218384, 7, 4255304), ncol = 2, byrow = T)

    #If both segments exist proceed normally
    if(noA == F & noB == F) {
        
        genomicA[,1] <- snpdata$POS[segmentcoordA[,1]] #Starting coords of A segments
        genomicA[,2] <- snpdata$POS[segmentcoordA[,2]] #End coords of A segments


        genomicB[,1] <- snpdata$POS[segmentcoordB[,1]] #Starting coords of B segments
        genomicB[,2] <- snpdata$POS[segmentcoordB[,2]] #End coords of B segments
        
        ##Start of the chromosome
        if(segmentcoordA[1,1] < segmentcoordB[1,1]) { genomicA[1,1] <- 1 } else { genomicB[1,1] <- 1 }

        ##End of the chromosome
        if(segmentcoordA[nsegmentA,2] > segmentcoordB[nsegmentB,2]) { genomicA[nsegmentA,2] <- ends[chr,2] } else { genomicB[nsegmentB,2] <- ends[chr,2] }

    }

    #If segment A does not exist
    if(noA == T) {
        genomicA[1,] <- 0
        genomicB[,1] <- 1
        genomicB[,2] <- ends[chr,2] }

    #If segment B does not exist
    if(noB == T) {
        genomicB[1,] <- 0
        genomicA[,1] <- 1
        genomicA[,2] <- ends[chr,2] }
    
    return(list("A" = genomicA, "B" = genomicB))
}



##This function uses inferred genomic coordinates to infer SNP genotypes
infer.SNPs <- function(genomicA, genomicB, parentA, parentB, currentname, allsnps) {

    SNPs <- NULL
    n.segA <- nrow(genomicA)
    n.segB <- nrow(genomicB)

    ##Looping over all A segments
    nameA <- paste("X", parentA, sep = "") ##Assuming that names are given as 2489 without the "X"
    indexA <- which(colnames(allsnps) == nameA) ##Col index of parent A
    for(k in 1:n.segA) {
        select.snp <- allsnps$pos >= genomicA[k,1] & allsnps$pos <= genomicA[k,2]
        temp <- allsnps[select.snp, c(4, indexA)]
        colnames(temp) <- c("pos", paste("X", currentname, sep = ""))
        SNPs <- rbind(SNPs, temp)
    }

    ##Looping over B segments
    nameB <- paste("X", parentB, sep = "") ##Assuming that names are given as 2489 without the "X"
    indexB <- which(colnames(allsnps) == nameB) ##Col index of parent B
    for(k in 1:n.segB) {
        select.snp <- allsnps$pos >= genomicB[k,1] & allsnps$pos <= genomicB[k,2]
        temp <- allsnps[select.snp, c(4, indexB)]
        colnames(temp) <- c("pos", paste("X", currentname, sep = ""))
        SNPs <- rbind(SNPs, temp)
    }


    #Now include SNPs that fall within area of unknown genotypes due to recombination breakpoints
    #And set these as missing data ("N")
    missing.snps <- !(allsnps$pos %in% SNPs$pos)
    snps.in.gaps <- allsnps[missing.snps, c(4, indexA, indexB)]
    ##If parents are the same, infer this for offspring, if they are different offspring data is missing
    inferred.snps <- ifelse(snps.in.gaps[,2] == snps.in.gaps[,3], snps.in.gaps[,2], "N")
    gap.snps <- data.frame(snps.in.gaps[,1], inferred.snps)
    #snps.in.gaps[,2] <- rep("N", nrow(snps.in.gaps)) #Set genotypes as missing data
    colnames(gap.snps) <- c("pos", paste("X", currentname, sep = ""))
    SNPs <- rbind(SNPs, gap.snps) #Combine all
    
    SNPs <- arrange(SNPs, pos) #Sort according to position

    return(SNPs)
}


#This function allows to easily check the parental genotypes of different individuals by chromosome
view.pgenot <- function(snpdata, genotypes, ind, chr) {

    filter <- snpdata$X.CHROM == paste("Supercontig_12.", chr, sep = "")
    snp.filt <- snpdata[filter,]
    df <- data.frame(chr = snp.filt$X.CHROM, position = snp.filt$POS, genotype = genotypes[[ind]]$GENOP[filter])
    df
}
