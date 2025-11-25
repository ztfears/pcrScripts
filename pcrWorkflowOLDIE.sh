#!/bin/bash --login

#SBATCH --time=23:59:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32G   # memory per CPU core
#SBATCH -J "pcrWorkflow"   # job name
###SBATCH --qos=bep8

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# Usage: sbatch pcrWorkflow.sh dirThatContainsStrainsFolder
# input file structure should be ../name/strains

CODE_DIR=$1
SCRIPTS_DIR=/grphome/grp_Assembly/pcrScripts
TIDY_ENV=/grphome/grp_Assembly/tidy


#Prokka
echo "Running submit_prokka.sh"
$SCRIPTS_DIR/submit_prokkaWF.sh $CODE_DIR/strains $CODE_DIR/prokka
#Wait for prokka jobs to finish
echo "Prokka finished!"

#Roary
echo "Running roary.sh"
sbatch --wait $SCRIPTS_DIR/roarySlurm.sh $CODE_DIR/prokka $CODE_DIR/roary
#Wait for roary job to finish
echo "Roary finished!"

#Find target genes
echo "Running findTargetGenesWrapper.sh"
$SCRIPTS_DIR/findTargetGenesWrapper.sh $CODE_DIR/roary/gene_presence_absence.csv $CODE_DIR/targetGenes
echo "findTargetGenes Finished!"

#Gather alleles
echo "Running gatherAllelesWrapper.sh"
$SCRIPTS_DIR/gatherAllelesWrapper.sh $CODE_DIR/targetGenes $CODE_DIR/prokka $CODE_DIR/alleles
echo "gatherAlleles Finished!"

#MAFFT
echo "Running mafft"
$SCRIPTS_DIR/submit_mafftWF.sh $CODE_DIR/alleles $CODE_DIR/mafft
echo "mafft Finished!"


#Either run entropyWrapper.sh or megacats
#if file starts with spec_, do:
#Run entropyWrapper.sh (Still need to test it works)
#mkdir -p $CODE_DIR/entropy

#if file starts with all_, do
#Run megaCATS
#mkdir -p $CODE_DIR/megaCATS
