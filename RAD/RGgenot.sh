#!/bin/bash

#Path to the reference genome
reference=~/Genomics/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta

#Paths to different analysis programs
path_to_picard_tools=~/picard-tools-1.119
path_to_GATK=~/GATK

#List of sample names (For family A)
#samplenames=~/Genomics/Neurospora/RAD/info/samplenamesfamA.txt
#samplenames=~/Genomics/Neurospora/RAD/info/samplenamesfamB.txt
#samplenames=~/Genomics/Neurospora/RAD/info/samplenamesfamC.txt
#samplenames=~/Genomics/Neurospora/RAD/info/samplenamesfamD.txt
#samplenames=~/Genomics/Neurospora/RAD/info/samplenamesfamE.txt
samplenames=~/Genomics/Neurospora/RAD/info/samplenamesfamG.txt

#First need to add read groups

while IFS= read -r sample; do

bamfile=~/Genomics/Neurospora/RAD/alignments/$sample.bam
outfile=~/Genomics/Neurospora/RAD/alignments/$sample.RG.bam

### Add or replace read-groups ###
#Note that sample names is used here in output and RGID!
java -Xmx4g -Djava.io.tmpdir='pwd'/tmp -jar $path_to_picard_tools/AddOrReplaceReadGroups.jar I=$bamfile O=$outfile SORT_ORDER=coordinate RGID=$sample RGLB=my_library RGPL=illumina RGSM=$sample RGPU=who_cares CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT TMP_DIR='pwd'/tmp

done < "$samplenames"

exit 0
