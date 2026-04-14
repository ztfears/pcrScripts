#!/bin/bash

### This script assumes your input CSV is just one column with all the unfiltered, unprocessed GenBank Accession numbers
### Output is a cleaned CSV that conforms to the correct standards and a .fna for each entry in the cleaned CSV


# Usage: First open a tmux session using: tmux new -s <session_name>
# Then run this: ./beaver.sh <path_to_inputCSV>
# Detach from the tmux session: Ctrl+B then D
# See what tmux sessions are available: tmux ls
# (if there are none when you first login, check the other login nodes by running ssh login01 or ssh login02 or ssh login03 or ssh login04, depending on which ones you've checked already)
# Check on tmux session: tmux attach -t <session_name>


#set -e
echo "Loading Conda..."
module load miniconda3

inCSV=$1
inDir=$(dirname $inCSV)

conda run -p /grphome/grp_Assembly/beaver Rscript --vanilla /grphome/grp_Assembly/pcrScripts/beaver/beaverCleaner.R $inCSV ##generates a temp.csv file in the inDir

for i in `cat $inDir/temp.csv`; do
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${i}&rettype=fasta&retmode=txt">$inDir/$i.fna;
done

### Remove empty files
find $inDir -type f ! -exec grep -q '[^[:space:]]' {} \; -delete
