# Neurospora nested association mapping population
This page describes a _Neurospora crassa_ nested association mapping population we have developed. We plan to add more strains to this population to keep improving it.

Last updated on 06.05.2024

Version 1.0

## Strains

This association mapping population contains strains that have been originally sampled from nature. Primarily from southeastern United States and the Caribbean basin. The natural strains have been obtained from Fungal Genetics Stock Center. In addition, we have created families, where we crossed some of the strains and the offspring of those crosses are also included in the NAM population.

Current stains numbers are:

|Sub-population|      n|    Parents|
|--------------|-------|-----------|
|Natural strains|    118|  |
|Family A|           94| 10948 × 10886    |
|Family B|           50| 10932 × 1165   |
|Family C|           50| 4498 × 8816   |
|Family D|           52| 3223 × 8845    |
|Family G|           70| 10904 × 851   |
|Total|              434|   |

## Genotyping

### Re-sequencing of the natural strains

The natural strains have had their genomes re-sequenced with Illumina, and either by us, or other groups. See [TPC GWAS preprint](https://www.biorxiv.org/content/10.1101/2024.04.29.591604v1) for details.

### RAD-sequencing the offspring

Since the genomes of the parents are known from sequencing, we genotyped the offspring using RAD-sequencing and then identified where recombination breakpoints occurred based on the RAD-genotypes. Then we could infer their full SNP genotypes from the SNP genotypes of the parents. See [TPC GWAS preprint](https://www.biorxiv.org/content/10.1101/2024.04.29.591604v1) for details.

### Genotype data file

The current genotype data file can be found at: [10.5281/zenodo.11120317](https://zenodo.org/records/11120317)

## Scripts
Like our previous work, we used a national supercluster (CSC) to process the short read sequencing data. The script files have been mainly written so that they work on the cluster. If you want to use them, you have to modify them so that they work in your environment. Nevertheless you can extract the GATK, BWA etc. commands and run them on your system.

The script files are in folder /resequencing

1. Map reads to to the reference genome using BWA: file /resequencing/mapping_array2.sh
2. Use haplotypecaller in the GATK pipeline to call genotypes for each sample (GVCF file)
3. Consolidate the samples into a database:
4. Jointly call genotypes from the sample database
5. Postprocess the resulting VCF files

