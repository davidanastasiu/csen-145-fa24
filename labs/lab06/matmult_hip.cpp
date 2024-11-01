#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdexcept>
#include <hip/hip_runtime.h>
#include "matmult.h"

#define BLOCK_SIZE 16

inline void gpuAssert(hipError_t code, const char *file, int line, bool abort=true)
{
   if (code != hipSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", hipGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#define gpuErrorCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }


/**
 * Matrix multiplication kernel
 * @param a GPU device pointer to an m X n matrix (A)
 * @param b GPU device pointer to a n X k matrix (B)
 * @param c GPU device output pointer to an m X k matrix (C) 
 * @param m nrows of A and C
 * @param n ncols of A and nrows of B
 * @param k ncols of B and C
 * 
 * Note:
 *     grid and block should be configured as:
 *         dim3 dimGrid((k + BLOCK_SIZE - 1) / BLOCK_SIZE, (m + BLOCK_SIZE - 1) / BLOCK_SIZE);
 *         dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
 */
__global__ void gpu_matmult(double *a, double *b, double *c, unsigned int m, unsigned int n, unsigned int k)
{ 
  // Implement matrix multiplication using HIP
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
  
  size_t memsize_A = sizeof(double) * nrows * ncols;
  size_t memsize_B = sizeof(double) * ncols * ncols2;
  size_t memsize_C = sizeof(double) * nrows * ncols2;

  // create host matrices
  auto A = create_mat(nrows, ncols, true);
  auto B = create_mat(ncols, ncols2, true);
  auto C = create_mat(nrows, ncols2, false);

  // create device matrices
  double *dA, *dB, *dC;
  gpuErrorCheck(hipMalloc((void **) &dA, memsize_A))
  gpuErrorCheck(hipMalloc((void **) &dB, memsize_B))
  gpuErrorCheck(hipMalloc((void **) &dC, memsize_C))

  // copy host memory to device
  gpuErrorCheck(hipMemcpy(dA, A, memsize_A, hipMemcpyHostToDevice))
  gpuErrorCheck(hipMemcpy(dB, B, memsize_B, hipMemcpyHostToDevice))

  // set up kernel execution parameters
  unsigned int grid_rows = ??; // calculate grid_rows
  unsigned int grid_cols = ??; // calculate grid_cols
  dim3 dimGrid(grid_cols, grid_rows);
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

  // execute kernel
  gpu_matmult <<< dimGrid, dimBlock >>> (dA, dB, dC, nrows, ncols, ncols2);

  // transfer results from device to host
  gpuErrorCheck(hipMemcpy(C, dC, memsize_C, hipMemcpyDeviceToHost))

  // optionally verify
  if(verify){
    // execute serial program on CPU
    auto C2 = create_mat(nrows, ncols2, false);
    serial_matmult(A, B, C2, nrows, ncols, ncols2);
    verify_result(C, C2, nrows, ncols2);
    free(C2);
  }

  // free memory
  gpuErrorCheck(hipFree(dA))
  gpuErrorCheck(hipFree(dB))
  gpuErrorCheck(hipFree(dC))
  free(A);
  free(B);
  free(C);

  return 0;
}
