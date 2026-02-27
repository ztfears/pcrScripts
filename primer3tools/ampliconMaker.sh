#!/bin/bash

# Usage: ampliconMaker.sh afaFile outputDir threshold

#SBATCH --time=00:30:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16G   # memory per CPU core
#SBATCH -J "ampliconMaker"   # job name

set -e

export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

input_dir=$1
output_dir=$2
threshold=$3

input_file=$(ls -1 $input_dir/*.afa | head -${SLURM_ARRAY_TASK_ID} | tail -1)

settingsFile=/grphome/grp_Assembly/pcrScripts/primer3_settings.txt

echo "Loading Conda ..."
module load miniconda3

stem=$(dirname "$input_file")
filenamer=$(basename "$input_file")
temp="${filenamer#*_}"
geneId="${temp%_*}"



echo "Processing $input_file..."
conda run --no-capture-output -p /grphome/grp_Assembly/seqinr bash -c "
  set -e
  Rscript --vanilla /grphome/grp_Assembly/pcrScripts/formEntWgtConsensus.R ${input_file} ${threshold}
  cat ${settingsFile} ${stem}/${geneId}*.txt | primer3_core > ${stem}/${geneId}_output.txt
"

mv $stem/$geneId* $output_dir

echo "Success!"
