#!/bin/bash --login
#submit with './mafft.sh <inputFile> <outputDir>'

set -e

unset MAFFT_BINARIES

file_name=$(basename "$1")

mafft --retree 1 --maxiterate 0 --memsave --thread 8 $1 > $2/$file_name.afa
