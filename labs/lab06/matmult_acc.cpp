#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdexcept>
#include <openacc.h>
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
void acc_matmult(double *a, double *b, double *c, unsigned int m, unsigned int n, unsigned int k) {
  // Implement matrix multiplication using OpenACC
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

  acc_device_t dev_type = acc_get_device_type();
  int ndevices = acc_get_num_devices(dev_type);
  cout << "Number of available devices: " << ndevices << endl;

  // create host matrices
  auto A = create_mat(nrows, ncols, true);
  auto B = create_mat(ncols, ncols2, true);
  auto C = create_mat(nrows, ncols2, false);

  acc_matmult(A, B, C, nrows, ncols, ncols2);

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