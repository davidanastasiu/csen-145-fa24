#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <algorithm>
#include <iostream>

void fillupRandomly(int *m, int size, unsigned int seed) {
    srand(seed);  
    for (int i = 0; i < size; i++) {
        m[i] = rand() % 10000;  
    }
}

void printArray(int *a, int size)
{
    for (int i = 0; i < size; i++)
        printf("%d ", a[i]);
    printf("\n");
}

int isSorted(int *a, int size)
{
    for (int i = 0; i < size - 1; i++)
        if (a[i] > a[i + 1])
            return 0;
    return 1;
}

int partition(int *a, int p, int r)
{
    int lt[r - p];
    int gt[r - p];
    int i;
    int j;
    int key = a[r];
    int lt_n = 0;
    int gt_n = 0;

    for (i = p; i < r; i++)
    {
        if (a[i] < a[r])
        {
            lt[lt_n++] = a[i];
        }
        else
        {
            gt[gt_n++] = a[i];
        }
    }

    for (i = 0; i < lt_n; i++)
    {
        a[p + i] = lt[i];
    }

    a[p + lt_n] = key;

    for (j = 0; j < gt_n; j++)
    {
        a[p + lt_n + j + 1] = gt[j];
    }

    return p + lt_n;
}

void quicksort(int *a, int p, int r)
{
    // Your Implementation
    // Your Implementation
    // Your Implementation
    // Your Implementation
    // Hint: OpenMP directives for parallelization
    
}

void merge(int *a, int size, int *b, int b_size, int *out)
{
    int i = 0, j = 0, k = 0;
    while (i < size && j < b_size)
    {
        if (a[i] <= b[j])
            out[k++] = a[i++];
        else
            out[k++] = b[j++];
    }
    while (i < size)
        out[k++] = a[i++];
    while (j < b_size)
        out[k++] = b[j++];
}

void gather_and_merge(int *local_data, int local_size, int N, int rank, int size)
{
    if (rank == 0)
    {
        int *sorted_data = (int *)malloc(N * sizeof(int));
        memcpy(sorted_data, local_data, local_size * sizeof(int));
        int offset = local_size;

        for (int i = 1; i < size; i++)
        {
            // Your Implementation
            // Your Implementation
            // Your Implementation
            // Your Implementation
            free(sorted_data);
            free(recv_data);
            sorted_data = temp;
            offset += recv_size;
        }

        assert(isSorted(sorted_data, N));

        printf("Array is sorted.\n");
        printf("Actual size of sorted data: %d\n", offset);

        free(sorted_data);
    }
    else
    {
        MPI_Send(&local_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(local_data, local_size, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = 1000000;
    int numThreads = (argc > 1) ? atoi(argv[1]) : 2;

    omp_set_num_threads(numThreads);

    int chunk_size = (N + size - 1) / size;
    int *local_data = (int *)malloc(chunk_size * sizeof(int));

    if (rank == 0)
    {
       // Your Implementation
      // Your Implementation
     // Your Implementation
    // Your Implementation
    }
    else
    {
        // Your Implementation
        // Your Implementation
        
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    
   
    // Your Implementation          
    // Your Implementation
    // Hint: OpenMP directives for parallelization
    quicksort(local_data, 0, chunk_size - 1);


    gather_and_merge(local_data, chunk_size, N, rank, size);

   //MPI_Barrier(MPI_COMM_WORLD); 
    double end_time = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "Total Sorting Time: " << (end_time - start_time) << " seconds" << std::endl;
    }

    free(local_data);
    MPI_Finalize();
    return 0;
}  
