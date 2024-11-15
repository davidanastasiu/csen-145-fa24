#!/bin/bash
#SBATCH --partition=mem          # Specify the partition
#SBATCH --nodes=2                # Number of nodes
#SBATCH --ntasks=16              # Number of tasks
#SBATCH --ntasks-per-node=8      # Number of tasks per node
#SBATCH --cpus-per-task=1        # Number of CPU cores per task
#SBATCH --mem=8G                 # Memory allocation
#SBATCH --time=01:00:00          # Max runtime (1 hour)            
#SBATCH --output=mpi_distribute_output.log # Log standard output to file
#SBATCH --error=mpi_distribute_error.log   # Log error output to file

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # Send email on job start, end, or failure
#SBATCH --mail-user=<your_email>@scu.edu  # Replace with your email address

# Load the OpenMPI module
module load OpenMPI/4.1.5
# Run the program
mpirun -np 16 ./lab08 /WAVE/projects/CSEN-145-Fa24/david/yelp/yelp.train.clu

