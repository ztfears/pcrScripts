#!/bin/bash

# Usage: entropyWrapper.sh afaFileDir outputDir

input_dir=$1
output_dir=$2

echo "Loading Conda"
module load miniconda3

for afa in $input_dir/*afa; do
	echo "Running entropyWF.R on $afa"
	conda run -p /grphome/grp_Assembly/seqinr Rscript --vanilla \
		/grphome/grp_Assembly/pcrScripts/entropyWF.R $afa
done;

mkdir -p $output_dir
mv $input_dir/*.csv $output_dir

echo "Success!" > $output_dir/success.txt
