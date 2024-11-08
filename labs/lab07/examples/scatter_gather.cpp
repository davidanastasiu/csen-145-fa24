#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    int root = 0;
    int data[8] = {10, 20, 30, 40,50,60,70,80};  // Data at root
    int recv_data;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    MPI_Scatter(data, 1, MPI_INT, &recv_data, 1, MPI_INT, root, MPI_COMM_WORLD);
    
    printf("Process %d received data: %d\n", rank, recv_data);
    
    recv_data *= 2;  // Modify data
    
    int gathered_data[4];
    MPI_Gather(&recv_data, 1, MPI_INT, gathered_data, 1, MPI_INT, root, MPI_COMM_WORLD);
    
    if (rank == root) {
        printf("Gathered data at root: ");
        for (int i = 0; i < size; i++) {
            printf("%d ", gathered_data[i]);
        }
        printf("\n");
    }
    
    MPI_Finalize();
    return 0;
}
