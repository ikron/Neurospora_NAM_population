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

outputdir=/scratch/project_2000350/genomics/natpop/genotyping/genoDB
inputdir=/scratch/project_2000350/genomics/natpop/genotyping
path_intervals=/scratch/project_2000350/genomics/Neurospora_reference/all.list

#Paths to programs
path_GATK=/projappl/project_2000350/Genomics/gatk-4.2.0.0


$path_GATK/gatk --java-options "-Xmx12g -Xms12g" GenomicsDBImport \
      -V $inputdir/10948.g.vcf.gz \
      -V $inputdir/10886.g.vcf.gz \
      -V $inputdir/10932.g.vcf.gz \
      -V $inputdir/1165.g.vcf.gz \
      -V $inputdir/4498.g.vcf.gz \
      -V $inputdir/8816.g.vcf.gz \
      -V $inputdir/3223.g.vcf.gz \
      -V $inputdir/8845.g.vcf.gz \
      -V $inputdir/10908.g.vcf.gz \
      -V $inputdir/847.g.vcf.gz \
      -V $inputdir/10904.g.vcf.gz \
      -V $inputdir/851.g.vcf.gz \
      -V $inputdir/1131.g.vcf.gz \
      -V $inputdir/8850.g.vcf.gz \
      -V $inputdir/8819.g.vcf.gz \
      -V $inputdir/4708.g.vcf.gz \
      -V $inputdir/4712.g.vcf.gz \
      -V $inputdir/6203.g.vcf.gz \
      -V $inputdir/4824.g.vcf.gz \
      -V $inputdir/8783.g.vcf.gz \
      -V $inputdir/8790.g.vcf.gz \
      -V $inputdir/3975.g.vcf.gz \
      -V $inputdir/10928.g.vcf.gz \
      -V $inputdir/10912.g.vcf.gz \
      -V $inputdir/4494.g.vcf.gz \
      -V $inputdir/3210.g.vcf.gz \
      -V $inputdir/10923.g.vcf.gz \
      -V $inputdir/10950.g.vcf.gz \
      -V $inputdir/10951.g.vcf.gz \
      -V $inputdir/10946.g.vcf.gz \
      -V $inputdir/3211.g.vcf.gz \
      -V $inputdir/10906.g.vcf.gz \
      -V $inputdir/P4452.g.vcf.gz \
      -V $inputdir/P4463.g.vcf.gz \
      -V $inputdir/P4468.g.vcf.gz \
      -V $inputdir/P4471.g.vcf.gz \
      -V $inputdir/P4476.g.vcf.gz \
      -V $inputdir/P4479.g.vcf.gz \
      -V $inputdir/10882.g.vcf.gz \
      -V $inputdir/10883.g.vcf.gz \
      -V $inputdir/10884.g.vcf.gz \
      -V $inputdir/10892.g.vcf.gz \
      -V $inputdir/10907.g.vcf.gz \
      -V $inputdir/10914.g.vcf.gz \
      -V $inputdir/10915.g.vcf.gz \
      -V $inputdir/10918.g.vcf.gz \
      -V $inputdir/10925.g.vcf.gz \
      -V $inputdir/10926.g.vcf.gz \
      -V $inputdir/10927.g.vcf.gz \
      -V $inputdir/10935.g.vcf.gz \
      -V $inputdir/10937.g.vcf.gz \
      -V $inputdir/10943.g.vcf.gz \
      -V $inputdir/10983.g.vcf.gz \
      -V $inputdir/3943.g.vcf.gz \
      -V $inputdir/P4489.g.vcf.gz \
      -L $path_intervals --genomicsdb-workspace-path $outputdir

