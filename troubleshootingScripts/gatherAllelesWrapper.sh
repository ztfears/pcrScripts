#!/bin/bash
#Usage: ./gatherAllelesWrapper.sh findTargetGenesOutputDir prokkaDir outputDir

set -e
ftgDir=$1
prokkaDir=$2
outDir=$3

mkdir -p $outDir

sampleCSV=$(ls -1 $ftgDir/*csv | head -n 1)

if [[ $sampleCSV == spec ]]; then
	fileTag="spec"
else
	fileTag=""
fi

tempFile=$ftgDir/"$fileTag"All.ffn

echo "Creating temp file"
#Combine all .ffn files for easy querying
cat $prokkaDir/"$fileTag"*.ffn > $tempFile

for csv in $ftgDir/*.csv; do /grphome/grp_Assembly/pcrScripts/gatherAlleles.sh $csv $tempFile; done

echo "Cleaning up temp file"
#Cleanup combined .ffn to save space over time
rm $tempFile

echo "Organizing utput"
mv $ftgDir/*.fna $outDir

#create success.txt
echo "Success!" > $outDir/success.txt
