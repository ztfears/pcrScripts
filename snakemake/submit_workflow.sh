#!/bin/bash
#submit with 'sbatch submit_workflow.sh <inputDir>'

#SBATCH --time=15:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8G   # memory per CPU core
#SBATCH -J "snakemake"   # job name
##SBATCH 

set -e

#INPUT_DIR=$1

INPUT_DIR=/home/zfears/Desktop/pcrProj/sm_eFaecium
SCRIPTS_DIR=/grphome/grp_Assembly/pcrScripts/snakemake/scripts

module load miniconda3
conda run -p /grphome/grp_Assembly/snakemake \
	snakemake --executor slurm --jobs 200 \
	--latency-wait 60 \
	--use-conda \
	--default-resources mem_mb_per_cpu=4000 \
	-s main.smk \
	--config input_dir="$INPUT_DIR" scripts_dir="$SCRIPTS_DIR"
