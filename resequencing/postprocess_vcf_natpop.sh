#!/bin/bash -l
#SBATCH --job-name=process_vcf
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

module load biokit
module load biojava/21
module load biopythontools

#Paths to programs
path_GATK=/projappl/project_2000350/Genomics/gatk-4.2.0.0
input=/scratch/project_2000350/genomics/natpop/genotyping

#Post processing of joint genotyping
#Merge vcf files
picard MergeVcfs I=$input/1.vcf.gz I=$input/2.vcf.gz I=$input/3.vcf.gz I=$input/4.vcf.gz I=$input/5.vcf.gz I=$input/6.vcf.gz I=$input/7.vcf.gz I=$input/8.vcf.gz I=$input/9.vcf.gz O=$input/all.samples.vcf.gz

#The combined .vcf is in file $all.samples.vcf.gz, the file needs to be indexed

#Index vcf file
$path_GATK/gatk --java-options "-Xmx4g" IndexFeatureFile -I $input/all.samples.vcf.gz
