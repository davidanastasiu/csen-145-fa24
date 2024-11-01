#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdexcept>
#include <omp.h>
#include "matmult.h"

/**
 * Matrix multiplication serial algorithm
 * 
 * @param a host device pointer to an m X n matrix (A)
 * @param b host device pointer to a n X k matrix (B)
 * @param c host device output pointer to an m X k matrix (C) 
 * @param m nrows of A and C
 * @param n ncols of A and nrows of B
 * @param k ncols of B and C
*/
void omp_gpu_matmult(double *a, double *b, double *c, unsigned int m, unsigned int n, unsigned int k) {
  #pragma omp target teams distribute parallel for simd map(to:a[0:m*n], b[0:n*k]) map(tofrom:c[0:m*k])
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < k; ++j)
    {
      double v = 0.0;
      for (int h = 0; h < n; ++h)
      {
        v += a[i * n + h] * b[h * k + j];
      }
      c[i * k + j] = v;
    }
  }
}

int main(int argc, char const *argv[])
{
    
  size_t nrows, ncols, ncols2;
  bool verify = argc > 4;
  
  if (argc < 4) {
    fprintf(stderr, "usage: matmult nrows ncols ncols2\n");
    return -1;
  }

  nrows = atoi(argv[1]);
  ncols = atoi(argv[2]);
  ncols2 = atoi(argv[3]);

  cout << "Matrix multiplication A(" << nrows << ", " << ncols << ") x B(" << ncols << ", " << ncols2 << ")" << endl;

  int nthreads = 1;
  #pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  cout << "CPU nthreads: " << nthreads << endl;

  int num_devices = omp_get_num_devices();
  cout << "Number of available devices: " << num_devices << endl;

  #pragma omp target 
  {
    if (omp_is_initial_device()) {
      printf("Running on host\n");    
    } else {
      int nteams= omp_get_num_teams(); 
      int nthreads= omp_get_num_threads();
      printf("Running on device with %d teams in total and %d threads in each team\n", nteams, nthreads);
    }
  }

  // create host matrices
  auto A = create_mat(nrows, ncols, true);
  auto B = create_mat(ncols, ncols2, true);
  auto C = create_mat(nrows, ncols2, false);

  omp_gpu_matmult(A, B, C, nrows, ncols, ncols2);

  // optionally verify
  if(verify){
    // execute serial program on CPU
    auto C2 = create_mat(nrows, ncols2, false);
    serial_matmult(A, B, C2, nrows, ncols, ncols2);
    verify_result(C, C2, nrows, ncols2);
    free(C2);
  }

  // free memory
  free(A);
  free(B);
  free(C);

  return 0;
}