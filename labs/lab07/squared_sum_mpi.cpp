#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define ARRAY_SIZE 32  // Total number of elements in the array

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int portion_size = ARRAY_SIZE / size;
    int* data = NULL;

    if (rank == 0) {
        // Root process generates the array
        data = (int*)malloc(ARRAY_SIZE * sizeof(int));
        for (int i = 0; i < ARRAY_SIZE; i++) {
            data[i] = i + 1;  // Example data: 1, 2, 3, ..., ARRAY_SIZE
        }
        printf("Original array:\n");
        for (int i = 0; i < ARRAY_SIZE; i++) {
            printf("%d ", data[i]);
        }
        printf("\n\n");
    }

    // Each process allocates space for its portion
    int* local_data = (int*)malloc(portion_size * sizeof(int));

    // Your Implementation (Hint: Use MPI )

   // Your Implementation

   // Your Implementation

   // Your Implementation

   // Your Implementation
    // Root process displays the result
    if (rank == 0) {
        printf("\nGlobal squared sum of array elements: %d\n", global_squared_sum);
    }

    // Cleanup
    if (rank == 0) free(data);
    free(local_data);

    MPI_Finalize();
    return 0;
}
