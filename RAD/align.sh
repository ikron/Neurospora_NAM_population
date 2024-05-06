#!/bin/bash

#Align RAD-tag reads to the reference

#file that contains sample names
#samplenames=~/Genomics/Neurospora/RAD/info/samplenames2.txt
samplenames=~/Genomics/Neurospora/RAD/info/samplenames.txt
reference=~/Genomics/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta

while IFS= read -r sample; do

    echo "Aligning sample" $sample
    fq_file=~/Genomics/Neurospora/RAD/cleaned/$sample.fq.gz
    bam_file=./alignments/$sample.temp.bam
    sorted_bam=./alignments/$sample

    #Aligning reads to reference using BWA
    bwa mem -M $reference $fq_file > tempfile.sam

    #Convert to .bam and sort
    samtools view -b -S tempfile.sam -o $bam_file

    samtools sort $bam_file $sorted_bam

    samtools index $sorted_bam.bam

    rm $bam_file

done < "$samplenames"

exit 0
