#!/bin/bash
#SBATCH --partition=hub          # Specify the partition
#SBATCH --nodes=1                # Number of nodes
#SBATCH --ntasks=1               # Number of tasks (1 task)
#SBATCH --cpus-per-task=8        # Number of CPU cores per task
#SBATCH --mem=2G                 # Memory allocation
#SBATCH --time=01:00:00          # Max runtime (1 hour)            
#SBATCH --output=cuda_output.log # Log standard output to file
#SBATCH --error=cuda_error.log   # Log error output to file

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # Send email on job start, end, or failure
#SBATCH --mail-user=your_email@example.com  # Replace with your email address

# Load the CUDA module
module load CUDA

# Compile the CUDA code
nvcc CUDA_problems.cu -o cuda_problems

# Run the compiled code
./cuda_problems

