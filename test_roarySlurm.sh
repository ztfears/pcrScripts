#!/bin/bash
###submit with 'sbatch roarySlurm.sh path_to_gff_files path_to_output_dir'

#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=20   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=64G   # memory per CPU core
#SBATCH -J "roary"   # job name
###SBATCH --qos=bep8

set -e
# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE


echo "Input GFF Files"
echo `$1/*.gff`



#start running the analysis
date
echo "Running roary..."
echo "roary -p 20 -e -n *.gff"

module load miniconda3
##running roary
conda run -p /nobackup/archive/grp/grp_Assembly/shared-conda-envs/roary \
        roary -e --mafft -p 20 -g 150000 -n -f $2 $1/*.gff

#Create success.txt
echo "Success!" > $2/success.txt

#call workflow
sbatch $SCRIPTS_DIR/betterWorkflow.sh $CODE_DIR
