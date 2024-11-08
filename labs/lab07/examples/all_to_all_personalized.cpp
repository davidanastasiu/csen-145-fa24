#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Each process prepares an array of data to send, where each element
    // is intended for a different process.
    int send_data[size];
    for (int i = 0; i < size; i++) {
        send_data[i] = rank * 10 + i;  // Unique data for each process
    }

    // Array to receive data from each process
    int recv_data[size];

    // Perform All-to-All communication
    MPI_Alltoall(send_data, 1, MPI_INT, recv_data, 1, MPI_INT, MPI_COMM_WORLD);

    // Display the results
    printf("Process %d received data: ", rank);
    for (int i = 0; i < size; i++) {
        printf("%d ", recv_data[i]);
    }
    printf("\n");

    MPI_Finalize();
    return 0;
}
