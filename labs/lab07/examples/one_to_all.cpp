#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank;
    int data = 100; // Data to broadcast
    int root = 0;   // Root process
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == root) {
        printf("Root process broadcasting data: %d\n", data);
    }
    
    MPI_Bcast(&data, 1, MPI_INT, root, MPI_COMM_WORLD);
    printf("Process %d received data: %d\n", rank, data);
    
    MPI_Finalize();
    return 0;
}
