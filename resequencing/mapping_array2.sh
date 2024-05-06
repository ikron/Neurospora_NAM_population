#!/bin/bash -l
#SBATCH --job-name=bwa_mapping
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-61
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

module load biokit

#Mapping all reads to the reference

#Path to the reference genome
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta
savedir=/scratch/project_2000350/genomics/natpop/alignments #Change outputdir when necessary

#List of sample names
samples=(1132 2229 3200 3212 3968 4713 4715 5914 7833 8784 8787 8789 8829 8848 8851 10881 10883 10884 10885 10887 10888 10889 10890 10891 10893 10894 10895 10896 10897 10898 10899 10900 10901 10902 10903 10905 10909 10910 10911 10916 10917 10919 10920 10921 10922 10929 10930 10931 10934 10936 10938 10939 10941 10942 10954 10982 P4451 P4457 P4459 P4472 P4486 P4496) #Make sure to include all samples


sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

fwd_fastq=/scratch/project_2000350/genomics/natpop/raw/${sample}_1.fq.gz
rev_fastq=/scratch/project_2000350/genomics/natpop/raw/${sample}_2.fq.gz

### Begin mapping reads ###
echo Mapping reads of $sample againts reference genome
#Note that since this would generate a lot of intermediate files I'm overwriting files up until final step.

### Using BWA
bwa mem $reference $fwd_fastq $rev_fastq > $savedir/${sample}.PE.sam

