#!/bin/bash 
# 
#SBATCH --job-name=someapp 
#SBATCH --output=DenseMatMultOmp.log 
# 
#SBATCH --partition=cmp
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00 
#
# Load the software to run the program 
module load GCC
cd ${HOME}/csen-145-fa24/labs/lab02
for t in 1 2 4 8 16; do
    export OMP_NUM_THREADS=${t}
    ./DenseMatMultOmp
done
