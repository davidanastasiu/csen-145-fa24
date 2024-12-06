#!/bin/bash
#SBATCH --partition=cmp          # Specify the partition, gpu or hub
#SBATCH --nodes=1                # Number of nodes
#SBATCH --ntasks-per-node=1            
#SBATCH --ntasks=1	            #total number of tasks
#SBATCH --cpus-per-task=8        # Number of CPU cores per task
#SBATCH --mem=4GB                 # Memory allocation
#SBATCH --time=01:10:00          # Max runtime (1 hour)      
#SBATCH --output=openblas_matmul_output.log # Log standard output to file
#SBATCH --error=openblas_matmul_error.log   # Log error output to file
# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # Send email on job start, end, or failure
#SBATCH --mail-user=rkachhadiya@scu.edu  # Replace with your email address

# Compile the OpenBlas  code
g++ -o openblas_matmul openblas_matmul.cpp -L/WAVE/users2/unix/rkachhadiya/.conda/envs/openBlas/lib -I/WAVE/users2/unix/rkachhadiya/.conda/envs/openBlas/include -lopenblas
export LD_LIBRARY_PATH=/WAVE/users2/unix/your_username/.conda/envs/openBlas/lib:$LD_LIBRARY_PATH
./openblas_matmul
