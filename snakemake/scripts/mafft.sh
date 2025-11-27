#!/bin/bash --login
#submit with 'sbatch mafft.sh <inputDir> <outputDir>'

#SBATCH --time=08:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16G   # memory per CPU core
#SBATCH -J "mafft"   # job name
##SBATCH 

set -e
# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

unset MAFFT_BINARIES


input_dir="$1"
file_name=$(ls -1 $input_dir | head -${SLURM_ARRAY_TASK_ID} | tail -1)
if [[ $file_name != *.fna ]]; then
        exit
fi

module load miniconda3
conda run -p /nobackup/archive/grp/grp_Assembly/shared-conda-envs/mafft \
	mafft --retree 1 --maxiterate 0 --memsave --thread 8 $1/$file_name > $2/$file_name.afa
