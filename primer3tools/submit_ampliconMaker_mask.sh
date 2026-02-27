#!/bin/bash
#Usage: ./submit_ampliconMaker_mask.sh inputDir labelsDir outputDir

input_dir=$1
labels_dir=$2
output_dir=$3
threshold=1.0

num_file=$(ls -1 $input_dir/*.afa | wc -l)
mkdir -p $output_dir
mkdir -p $output_dir/slurm_results

sbatch --array 1-$num_file -o $output_dir/slurm_results/slurm-%A_%a.out /grphome/grp_Assembly/pcrScripts/ampliconMaker_mask.sh $input_dir $labels_dir $output_dir $threshold
