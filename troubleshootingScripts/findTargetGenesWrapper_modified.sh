#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=4G   # memory per CPU core
#SBATCH -J "ftg_mod"



# Usage: ./findTargetGenesWrapper_modified.sh /path/to/gene_presence_absence.csv outputDir

set -e
echo "Loading Conda"
module load miniconda3

echo "Running findTargetGenes_modified.R"

mkdir -p $2
conda run -p /grphome/grp_Assembly/tidy Rscript --vanilla /grphome/grp_Assembly/pcrScripts/troubleshootingScripts/findTargetGenes_modified.R $1 $2


#create success.txt
echo "Success!" > $2/success.txt
