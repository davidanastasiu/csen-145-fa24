#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>  // For std::sort

#define ARRAY_SIZE 32  // Total number of elements to be sorted
#define ROOT 0

void generate_data(int* array, int size) {
    srand(time(NULL));
    for (int i = 0; i < size; i++) {
        array[i] = rand() % 100; // Random values between 0 and 99
    }
}

void merge(int* result, int* left, int left_size, int* right, int right_size) {
    int i = 0, j = 0, k = 0;
    while (i < left_size && j < right_size) {
        if (left[i] < right[j]) {
            result[k++] = left[i++];
        } else {
            result[k++] = right[j++];
        }
    }
    while (i < left_size) result[k++] = left[i++];
    while (j < right_size) result[k++] = right[j++];
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int* global_data = NULL;
    int chunk_size = ARRAY_SIZE / size;

    if (rank == ROOT) {
        global_data = (int*)malloc(ARRAY_SIZE * sizeof(int));
        generate_data(global_data, ARRAY_SIZE);

        printf("Original unsorted array:\n");
        for (int i = 0; i < ARRAY_SIZE; i++) {
            printf("%d ", global_data[i]);
        }
        printf("\n\n");
    }

    // Each process allocates memory for its chunk
    int* local_data = (int*)malloc(chunk_size * sizeof(int));

   // Your Implementation (Hint: Use MPI and use std::sort )

   // Your Implementation

   // Your Implementation

   // Your Implementation

   // Your Implementation

    // Root process merges the sorted chunks
    if (rank == ROOT) {
        // Iteratively merge pairs of sorted chunks
        int* temp_array = (int*)malloc(ARRAY_SIZE * sizeof(int));
        int current_size = chunk_size;

        while (current_size < ARRAY_SIZE) {
            for (int i = 0; i < ARRAY_SIZE; i += 2 * current_size) {
                int left_size = current_size;
                int right_size = (i + 2 * current_size <= ARRAY_SIZE) ? current_size : (ARRAY_SIZE - i - current_size);
                
                if (right_size > 0) {
                    merge(temp_array + i, sorted_data + i, left_size, sorted_data + i + left_size, right_size);
                } else {
                    std::copy(sorted_data + i, sorted_data + i + left_size, temp_array + i);
                }
            }
            std::swap(sorted_data, temp_array);
            current_size *= 2;
        }
        free(temp_array);

        printf("Globally sorted array:\n");
        for (int i = 0; i < ARRAY_SIZE; i++) {
            printf("%d ", sorted_data[i]);
        }
        printf("\n");
        
        free(global_data);
        free(sorted_data);
    }

    free(local_data);
    MPI_Finalize();
    return 0;
}
