#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size, local_data, global_sum;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    local_data = rank + 1;  // Each processor has different data

    // Print the local data at each process
    printf("Process %d out of %d has local data: %d\n", rank, size, local_data);
    
    // Perform the reduction operation to compute the global sum at the root process
    MPI_Reduce(&local_data, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    // Only the root process prints the final result
    if (rank == 0) {
        printf("Total sum of data across all processes: %d\n", global_sum);
    }
    
    MPI_Finalize();
    return 0;
}
