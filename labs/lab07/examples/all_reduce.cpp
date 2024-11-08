#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Each process starts with a unique integer value, here simply rank + 1
    int local_value = rank + 1;

    // Variables to store the results of the all-reduce operations
    int global_sum, global_max;

    // Perform Allreduce to compute the global sum of local values
    MPI_Allreduce(&local_value, &global_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // Perform Allreduce to compute the global maximum of local values
    MPI_Allreduce(&local_value, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // Display the results
    printf("Process %d: local value = %d, global sum = %d, global max = %d\n", rank, local_value, global_sum, global_max);

    MPI_Finalize();
    return 0;
}
