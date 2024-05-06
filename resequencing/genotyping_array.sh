#!/bin/bash -l
#SBATCH --job-name=haplotypecaller
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

module load biokit
module load biopythontools

#Using haplotype caller to call genotypes one by one

#Path to the reference genome
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta
inputdir=/scratch/project_2000350/genomics/natpop/alignments
outputdir=/scratch/project_2000350/genomics/natpop/genotyping

#Paths to programs
path_GATK=/projappl/project_2000350/Genomics/gatk-4.2.0.0

samples=(1132 2229 3200 3212 3968 4713 4715 5914 7833 8784 8787 8789 8829 8848 8851 10881 10883 10884 10885 10887 10888 10889 10890 10891 10893 10894 10895 10896 10897 10898 10899 10900 10901 10902 10903 10905 10909 10910 10911 10916 10917 10919 10920 10921 10922 10929 10930 10931 10934 10936 10938 10939 10941 10942 10954 10982 P4451 P4457 P4459 P4472 P4486 P4496)


sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

$path_GATK/gatk --java-options "-Xmx8g" HaplotypeCaller  \
   -R $reference -I $inputdir/$sample.RG.bam --sample-ploidy 2 -mbq 10 --output-mode EMIT_ALL_CONFIDENT_SITES -O $outputdir/$sample.g.vcf.gz -ERC GVCF 
