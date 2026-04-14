#!/bin/bash

# Usage: gatherAlleles.sh CsvW/ProkkaID tempFile outDir

set -e
tempFile=$2
prokkaIDCsv="$1"

filename=$(basename "$prokkaIDCsv")
base="${filename%.csv}"
outDir=$3

#Ensure UNIX Csv formatting
dos2unix -q $prokkaIDCsv

echo "Loading seqtk"
module load seqtk

echo "Pull Genes by ID"
#Pull Genes by ID using seqtk
seqtk subseq $tempFile $prokkaIDCsv > $outDir/"$base"_alleles.fna

echo "Finished!"
