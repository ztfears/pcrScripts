#!/bin/bash
#Usage: ./submit_ampliconMaker.sh inputDir outputDir

input_dir=$1
output_dir=$2
threshold=1.0

num_file=$(ls -1 $input_dir/*.afa | wc -l)
mkdir -p $output_dir
mkdir -p $output_dir/slurm_results

sbatch --array 1-$num_file -o $output_dir/slurm_results/slurm-%A_%a.out /grphome/grp_Assembly/pcrScripts/ampliconMaker.sh $input_dir $output_dir $threshold
