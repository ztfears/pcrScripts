#!/bin/bash

# Usage: ./findTargetGenesWrapper.sh path_to_genePresenceAbsence.csv outputDir

set -e
echo "Loading Conda"
module load miniconda3

echo "Running findTargetGenes_modified.R"

mkdir -p $2
conda run -p /grphome/grp_Assembly/tidy Rscript --vanilla /grphome/grp_Assembly/pcrScripts/troubleshootingScripts/findTargetGenes.R $1 $2


#create success.txt
echo "Success!" > $2/success.txt
