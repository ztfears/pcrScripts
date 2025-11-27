#!/bin/bash
###submit with './prokka.sh inputfile output_dir'

set -e

#start running the analysis
date
echo "Running prokka..."
echo "prokka file..."

##running prokka
input_dir=$(dirname "$1")
file_name=$(basename "$1")
if [[ $file_name == *.csv ]]; then
	exit
fi

#Checking if file is empty
size=$(du $input_dir/$file_name | tr '\t' '\n'| head -1)
if [[ $size -le 1 ]]; then
	exit
fi


output_dir="$2"
echo $output_dir
echo $input_dir/$file_name

#module load miniconda3
#conda run -p /nobackup/archive/grp/grp_Assembly/shared-conda-pkgs/envs/prokka \
	prokka --outdir $output_dir --prefix $file_name --force $input_dir/$file_name

