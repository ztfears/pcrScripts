#!/bin/bash
#Usage: ./submit_mafft.sh inputDir outputDir

input_dir=$1
if [ -z "$input_dir" ]; then
	echo "Pls provide input directory you fetcher"
	exit
fi

output_dir=$2
if [ -z "$output_dir" ]; then
	echo "Pls provide output directory you fetcher"
	exit
fi

num_file=$(ls -1 $input_dir | wc -l)
mkdir -p $output_dir/slurm_results

sbatch --array 1-$num_file -o $output_dir/slurm_results/slurm-%A_%a.out /grphome/grp_Assembly/pcrScripts/troubleshootingScripts/mafft.sh $input_dir $output_dir
