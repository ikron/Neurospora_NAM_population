#!/bin/bash

#Enable devtoolset-3 first to get gcc 4.9
#scl enable devtoolset-3 bash

#Note! devtoolset-3 has to be enabled to run stacks. However, run the command above before executing any scripts because the scl enable foo bash runs a new bash instance

#Cleanup and demultiplex RAD-seq reads

#Plate 1
rawdir=~/Genomics/Neurospora/RAD/raw/plate1/

barcodes_file=~/Genomics/Neurospora/RAD/info/plate1.tsv

process_radtags -p $rawdir -b $barcodes_file -o ./ -e pstI --inline_null -c -q -r &> process_radtags.plate1.oe

#Plate 2
rawdir=~/Genomics/Neurospora/RAD/raw/plate2/
barcodes_file=~/Genomics/Neurospora/RAD/info/plate2.tsv

process_radtags -p $rawdir -b $barcodes_file -o ./ -e pstI --inline_null -c -q -r &> process_radtags.plate2.oe

#Plate 3
rawdir=~/Genomics/Neurospora/RAD/raw/plate3/
barcodes_file=~/Genomics/Neurospora/RAD/info/plate3.tsv

process_radtags -p $rawdir -b $barcodes_file -o ./ -e pstI --inline_null -c -q -r &> process_radtags.plate3.oe

#Plate 4
rawdir=~/Genomics/Neurospora/RAD/raw/plate4/
barcodes_file=~/Genomics/Neurospora/RAD/info/plate4.tsv

process_radtags -p $rawdir -b $barcodes_file -o ./ -e pstI --inline_null -c -q -r &> process_radtags.plate4.oe
