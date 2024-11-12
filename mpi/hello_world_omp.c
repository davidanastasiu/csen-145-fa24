#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int npes;
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    // Get the rank of the process
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Get the name of the process
    char pname[MPI_MAX_PROCESSOR_NAME];
    int sz;
    MPI_Get_processor_name(pname, &sz);

    // Get/set the number of threads
    if(argc > 1){
        int nthreads = atoi(argv[1]);
        omp_set_num_threads(nthreads);
    }

    // Print the hello world message
#pragma omp parallel
{
    int tid = omp_get_thread_num();
    printf("    From thread %d of process %s, rank %d, out of %d processes, Hello World!\n",
           tid, pname, myrank, npes);
}

    // Finalize the MPI environment.
    MPI_Finalize();
}
