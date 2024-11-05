#include "saxpy.h"
#include <openacc.h>

double saxpy(int n,
  float a,
  float *x,
  float *restrict y)
{
  double time;
  #pragma acc data copyin(x[0:N]) copy(y[0:N])
  {
  saxpy_timer timer;
  #pragma acc parallel loop
  for (int i = 0; i < n; ++i)
    y[i] = a*x[i] + y[i];
  time = timer.elapsed_msec();
  }
  return time;
}

int main() {
   float *x = new float[N], *y = new float[N];

   for (int i = 0; i < N; ++i) {
      x[i] = XVAL;
      y[i] = YVAL;
   }
   std::cout << "N: " << N << std::endl;

   saxpy_timer timer;
   auto elapsed2 = saxpy(N, AVAL, x, y);
   auto elapsed = timer.elapsed_msec();
   std::cout << "Elapsed: " << elapsed << " ms\n";
   std::cout << "Compute only: " << elapsed2 << " ms\n";

   saxpy_verify(y);
   delete[] x;
   delete[] y;
   return 0;
}

