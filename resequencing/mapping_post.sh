#!/bin/bash -l
#SBATCH --job-name=mapping_post
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-61
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

module load biokit

#Set temp directory
TMP_DIR=/scratch/project_2000350/genomics/

#Process the alignments and convert SAM files to BAM files and 

#Path to the reference genome
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta
savedir=/scratch/project_2000350/genomics/natpop/alignments #Change outputdir when necessary

#List of sample names
samples=(1132 2229 3200 3212 3968 4713 4715 5914 7833 8784 8787 8789 8829 8848 8851 10881 10883 10884 10885 10887 10888 10889 10890 10891 10893 10894 10895 10896 10897 10898 10899 10900 10901 10902 10903 10905 10909 10910 10911 10916 10917 10919 10920 10921 10922 10929 10930 10931 10934 10936 10938 10939 10941 10942 10954 10982 P4451 P4457 P4459 P4472 P4486 P4496)

sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

#Convert to a .bam file
samtools view -S -b $savedir/${sample}.PE.sam > $savedir/${sample}.PE.bam

#Sorting reads
samtools sort $savedir/${sample}.PE.bam -o $savedir/${sample}.PE.sorted.bam

#index bam
samtools index $savedir/${sample}.PE.sorted.bam

### Fix mate information ###
picard FixMateInformation VALIDATION_STRINGENCY=LENIENT INPUT=$savedir/${sample}.PE.sorted.bam OUTPUT=$savedir/${sample}.mate_fixed.bam

#Sort the resulting bam
samtools sort $savedir/${sample}.mate_fixed.bam -o $savedir/${sample}.mate_fixed.sorted.bam

### Add or replace read-groups ###
#Note that sample names is used here in output and RGID!
picard AddOrReplaceReadGroups I=$savedir/${sample}.mate_fixed.sorted.bam O=$savedir/$sample.RG.bam SORT_ORDER=coordinate RGID=$sample RGLB=my_library RGPL=illumina RGSM=$sample RGPU=who_cares CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

#Remove temp files so that they don't accumulate and take space
rm $savedir/$sample.PE.sam 
rm $savedir/$sample.PE.bam
rm $savedir/$sample.PE.sorted.bam
rm $savedir/$sample.mate_fixed.bam
rm $savedir/$sample.mate_fixed.sorted.bam

