#!/bin/bash -l
#SBATCH --job-name=dbimport
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=26G
#SBATCH --mail-type=END
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

module load biokit
module load biopythontools

outputdir=/scratch/project_2000350/genomics/natpop/genotyping/genoDB
inputdir=/scratch/project_2000350/genomics/natpop/genotyping
path_intervals=/scratch/project_2000350/genomics/Neurospora_reference/all.list

#Paths to programs
path_GATK=/projappl/project_2000350/Genomics/gatk-4.2.0.0

$path_GATK/gatk --java-options "-Xmx12g -Xms12g" GenomicsDBImport \
      -V $inputdir/5910.g.vcf.gz \
      -V $inputdir/4730.g.vcf.gz \
      -V $inputdir/1133.g.vcf.gz \
      -V $inputdir/4716.g.vcf.gz \
      -L $path_intervals --genomicsdb-update-workspace-path $outputdir
