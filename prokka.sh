#!/bin/bash --login
###submit with 'sbatch prokka.sh path_to_fasta_files'

#SBATCH --time=00:15:00   # walltime
#SBATCH --ntasks=10   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32G   # memory per CPU core
#SBATCH -J "prokka"   # job name
###SBATCH --qos=bep8

set -e
# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#initialize conda and activate environment
#source /fslhome/bep8/anaconda3/bin/activate #local
#module load group-conda

module load miniconda3
#source /nobackup/scratch/grp/grp_Assembly/shared-conda-pkgs/bin/activate

#start running the analysis
date
echo "Running prokka..."
echo "prokka file..."

##running prokka
input_dir="$1"
file_name=$(ls -1 $input_dir | head -${SLURM_ARRAY_TASK_ID} | tail -1)
if [[ $file_name == *.csv ]]; then
	exit
fi

#Checking if file is empty
size=$(du $1/$file_name | tr '\t' '\n'| head -1)
if [[ $size -le 1 ]]; then
	exit
fi


output_dir="$2"
echo $output_dir
echo $input_dir/$file_name
conda run -p /nobackup/archive/grp/grp_Assembly/shared-conda-pkgs/envs/prokka \
	prokka --outdir $output_dir --prefix $file_name --force $input_dir/$file_name

