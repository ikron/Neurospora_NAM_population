#!/bin/bash -l
#SBATCH --job-name=GVCF_array
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=48G
#SBATCH --array=1-9
#SBATCH --mail-type=END
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

module load biokit
module load biopythontools

#Using genotype GVCF to consolidate genotype calls in a database into one vcf file

#Path to the reference genome
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta
outputdir=/scratch/project_2000350/genomics/natpop/genotyping
variantdb=/scratch/project_2000350/genomics/natpop/genotyping/genoDB
intervals=/scratch/project_2000350/genomics/mutacc/genotyping/interval_list #List of genomic intervals, used in parallelization

#Paths to programs
path_GATK=/projappl/project_2000350/Genomics/gatk-4.2.0.0


#Note annoying syntax '-V gendb:///home/path/myDB' is correct three (3) forward slashes!
$path_GATK/gatk --java-options "-Xmx48g" GenotypeGVCFs  \
   -R $reference -V gendb://$variantdb --sample-ploidy 2 -all-sites -L $intervals/${SLURM_ARRAY_TASK_ID}.list --tmp-dir $outputdir -O $outputdir/${SLURM_ARRAY_TASK_ID}.vcf.gz

#Note that resulting vcf.files need to be combined and indexed

