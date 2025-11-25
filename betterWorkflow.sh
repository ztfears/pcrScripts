#!/bin/bash --login

#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32G   # memory per CPU core
#SBATCH -J "betterWorkflow"   # job name
###SBATCH --qos=bep8

set -e
# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# Usage: sbatch betterWorkflow.sh dirThatContainsStrainsFolder
# input file structure should be ../name/strains

export CODE_DIR=$1
export SCRIPTS_DIR=/grphome/grp_Assembly/pcrScripts

echo "$CODE_DIR"

#Prokka
if [ ! -f $CODE_DIR/prokka/success.txt ]; then
	echo "Running submit_prokka.sh"
	$SCRIPTS_DIR/submit_prokkaWF.sh $CODE_DIR/strains $CODE_DIR/prokka
fi

#Roary
if [ ! -f $CODE_DIR/roary/success.txt ]; then
	echo "Prokka finished!"
	echo "Running roary.sh"
	sbatch -o $CODE_DIR/roary/slurm_results/slurm-%j.out  $SCRIPTS_DIR/roarySlurm.sh $CODE_DIR/prokka $CODE_DIR/roary
	exit
fi
echo "Roary finished!"


#Find target genes
if [ ! -f $CODE_DIR/targetGenes/success.txt ]; then
	echo "Running findTargetGenesWrapper.sh"
	$SCRIPTS_DIR/findTargetGenesWrapper.sh $CODE_DIR/roary/gene_presence_absence.csv $CODE_DIR/targetGenes
	echo "findTargetGenes Finished!"
fi

#Gather alleles
if [ ! -f $CODE_DIR/alleles/success.txt ]; then
	echo "Running gatherAllelesWrapper.sh"
	$SCRIPTS_DIR/gatherAllelesWrapper.sh $CODE_DIR/targetGenes $CODE_DIR/prokka $CODE_DIR/alleles
	echo "gatherAlleles Finished!"
fi

#MAFFT
if [ ! -f $CODE_DIR/mafft/success.txt ]; then
	echo "Running mafft"
	$SCRIPTS_DIR/submit_mafftWF.sh $CODE_DIR/alleles $CODE_DIR/mafft
	echo "mafft Finished!"
fi

##If species-only genes were found, run entropy
##Entropy
#outFile=$(ls -1 /home/zfears/Desktop/pcrProj/sPneumo/mafft/ | head -2 | tail -1)
#echo $outFile
#echo `$outFile == spec_*`

#if [[ $outFile == spec_* ]]; then
#	if [ ! -f $CODE_DIR/entropy/success.txt ]; then
#		echo "Running entropy"
#		$SCRIPTS_DIR/entropyWrapper.sh $CODE_DIR/mafft 200 20 $CODE_DIR/entropy
#		echo "entropy Finished!"
#	fi
#fi

