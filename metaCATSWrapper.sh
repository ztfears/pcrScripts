#!/bin/bash

# Usage: metaCATSWrapper.sh afa&csvFileDir outputDir

#all_nusG_alleles.fna.afa
#all_nusG_locVstr.csv

input_dir=$1
output_dir=$2

echo "Loading Conda"
module load miniconda3

for afa in $input_dir/*afa; do
	echo "Running entropyWF.R on $afa"
	## strip afa of ending up to gene, then combine to make locVstr for thingey
	stem=${afa%_*}
	meta="${stem}_locVstr.csv"
	conda run -p /grphome/grp_Assembly/seqinr Rscript --vanilla \
		/grphome/grp_Assembly/pcrScripts/metaCATS-PCR.R $afa $meta
done;

mkdir -p $output_dir
mv $input_dir/*.csv $output_dir

echo "Success!" > $output_dir/success.txt
