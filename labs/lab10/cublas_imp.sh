#!/bin/bash
#SBATCH --partition=hub          # Specify the partition, gpu or hub
#SBATCH --nodes=1                # Number of nodes
#SBATCH --ntasks-per-node=1            
#SBATCH --ntasks=1	            #total number of tasks
#SBATCH --cpus-per-task=8        # Number of CPU cores per task
#SBATCH --mem=4GB                 # Memory allocation
#SBATCH --time=00:10:00          # Max runtime (1 hour)      
#SBATCH --output=cublas_matmul_output.log # Log standard output to file
#SBATCH --error=cublas_matmul_error.log   # Log error output to file
# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # Send email on job start, end, or failure
#SBATCH --mail-user=rkachhadiya@scu.edu  # Replace with your email address

module load CUDA
nvcc -o cublas_matmul cublas_matmul.cu -lcublas
./cublas_matmul
