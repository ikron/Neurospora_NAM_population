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

### Analysis of resequencing data

1. Map reads to to the reference genome using BWA: file /resequencing/mapping_array2.sh
2. Postprocess the mapping files: file /resequencing/mapping_post.sh
3. Use haplotypecaller in the GATK pipeline to call genotypes for each sample (GVCF file): file /resequencing/genotyping_array.sh
4. Consolidate the samples into a database: file/resequencing/db_import_natpop.sh and file /resequencing/db_import_update.sh
5. Jointly call genotypes from the sample database: file /resequencing/genotypeGVCF_natpop.sh
6. Postprocess the resulting VCF files: file /resequencing/postprocess_vcf_natpop.sh
7. Further VCF post processing

  -  Extract .gz file

  -  Remove INFO.LEAC, INFO:MLEAF columns from vcf file using bcftools

    ~/bcftools-1.11/bcftools annotate -x INFO/MLEAC,INFO/MLEAF all.samples.vcf -o all.samples.rm.vcf

  -  For nat pop also run (PID column too long, and Neurospora is haploid anyway)

    ~/bcftools-1.11/bcftools annotate -x FORMAT/PID all.samples.rm.vcf -o all.samples.rm2.vcf

    ~/bcftools-1.11/bcftools annotate -x INFO/MLEAC,INFO/MLEAF,FORMAT/PID all.samples.vcf -o all.samples.natpop_comp.vcf

 8. Then make an indexed database from the resulting VCF file using wormtable

  -  Convert to wormtable database

    vcf2wt --truncate all.samples.rm.vcf allsamples.wt

    wtadmin add allsamples.wt CHROM+POS

  -  Can look what are the contents of wormtable columns with wtadmin show allsamples.wt

9. Filter SNPs called by the GATK pipeline: file /resequencing/natpop_genotypes.py
   - The output file from this step will be used in inferring genotypes of the offspring based on RAD-sequencing

### Analysis of RAD-sequencing data

10. Clean and demultiplex RAD reads: file /RAD/cleanup.sh

Note that the .fastq files uploaded to SRA are the demultiplexed fastq files for each sample

11. Align RAD-reads to the reference genome: file /RAD/align.sh
   - This uses the samplenames.txt to loop over all samples

12. Add read groups to the samples: file /RAD/RGgenot.sh
   - This uses the samplenamesfamX.txt files to process the different families

13. Local indel realignment for all samples in the same family simultaneously: file /RAD/indelrealign.sh
   - This uses the famXbams.list files to group families

14. Genotyping family bams with GATK Unified Genotyper: file /RAD/GATKgenot.sh
    - This genotypes the family based bam files and produces a .vcf file for the family. These are used as input files in the next step

15. Process the vfamily-wise .vcf files and convert SNP calls into parent calls. Then infer the full SNP data for the offspring based on recombination breakpoints. File: /RAD/RADgeno.R
   - This first produces family-wise genotype calls, and then these are subsequently used together with the hapmap data file from step 9. to produce the final genotype file
