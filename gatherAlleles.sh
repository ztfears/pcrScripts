#!/bin/bash

# Usage: gatherAlleles.sh CsvW/ProkkaID tempFile

set -e
tempFile=$2
prokkaIDCsv="$1"
base="${prokkaIDCsv%.csv}"

#Ensure UNIX Csv formatting
dos2unix -q $prokkaIDCsv

echo "Loading seqtk"
module load seqtk

echo "Pull Genes by ID"
#Pull Genes by ID using seqtk
seqtk subseq $tempFile $prokkaIDCsv > "$base"_alleles.fna

echo "Finished!"
