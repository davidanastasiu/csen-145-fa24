#include <mpi.h>
#include <stdio.h>

#define ARRAY_SIZE 4  // Each process holds 4 elements

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != ARRAY_SIZE) {
        if (rank == 0) {
            printf("This example requires %d processes.\n", ARRAY_SIZE);
        }
        MPI_Finalize();
        return -1;
    }

    // Each process has a unique "row" of the matrix
    int send_data[ARRAY_SIZE];
    for (int i = 0; i < ARRAY_SIZE; i++) {
        send_data[i] = rank * 10 + i;  // Unique data for each process
    }

    // Define the receive count for each process (each process receives 1 element)
    int recv_counts[ARRAY_SIZE] = {1, 1, 1, 1};

    // Buffer to hold the result of MPI_Reduce_scatter
    int recv_data;

    // Perform Reduce-Scatter with a sum operation
    MPI_Reduce_scatter(send_data, &recv_data, recv_counts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // Display the results
    printf("Process %d received the sum of column %d: %d\n", rank, rank, recv_data);

    MPI_Finalize();
    return 0;
}
