#include <iostream>
#include <omp.h>
#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */
#include <cstring>  /* strcasecmp */

using namespace std;

int *random_vector(int size)
{
    auto v = new int[size];
    unsigned int seed = rand();
    for (int i = 0; i < size; ++i)
    {
        v[i] = rand_r(&seed);
    }
    return v;
}

void dotp(int n, int* v1, int* v2)
{
    /* find dot product */
    double dotp = 0.0;

    #pragma omp parallel for reduction(+ : dotp)
    for (int i = 0; i < n; ++i)
    {
        dotp += v1[i] * v2[i];
    }
    cout << "dotp = " << dotp << endl;
}


int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Invalid options." << endl
             << "<program> <num_elements> [-t <num_threads>]" << endl;
        exit(1);
    }
    int n = atoi(argv[1]);
    if (argc == 4 && strcasecmp(argv[2], "-t") == 0)
    {
        int nthreads = atoi(argv[3]);
        omp_set_num_threads(nthreads);
    }

    /* initialize random seed: */
    srand(time(NULL));
    
    auto vec1 = random_vector(n);
    auto vec2 = random_vector(n);

    #pragma omp parallel
    {
        int nt = omp_get_num_threads();
        #pragma omp single
        cout << "Executing with " << nt << " threads." << endl;
    }

    double start;
    double end;
    start = omp_get_wtime();
    dotp(n, vec1, vec2);
    end = omp_get_wtime();
    cout << "Work took " << end - start << " seconds." << endl;


    delete[] vec1;
    delete[] vec2;

    return 0;
}
