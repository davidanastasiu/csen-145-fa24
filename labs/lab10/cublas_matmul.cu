#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cmath>  


bool compareMatrices(float* a, float* b, int n, float tol = 1e-3) {
    for (int i = 0; i < n * n; i++) {
        if (std::fabs(a[i] - b[i]) > tol) {
            return false; 
        }
    }
    return true;  
}

// Error handling macro
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}


void matrixMultiplyCPU(float* a, float* b, float* c, int n) {
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            float sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += a[row * n + k] * b[k * n + col];
            }
            c[row * n + col] = sum;
        }
    }
}

int main() {
    int n = 1000;
    float *a, *b, *c, *c_serial;
    float *d_a, *d_b, *d_c;
    srand(42);
    // Allocate host memory
    a = (float*)malloc(n * n * sizeof(float));
    b = (float*)malloc(n * n * sizeof(float));
    c = (float*)malloc(n * n * sizeof(float));
    c_serial = (float*)malloc(n * n * sizeof(float));

   
    for(int i = 0; i < n * n; i++) {
        a[i] = static_cast<float>(rand()) / RAND_MAX;
        b[i] = static_cast<float>(rand()) / RAND_MAX;
    }

    
    cudaMalloc((void**)&d_a, n * n * sizeof(float));
    cudaMalloc((void**)&d_b, n * n * sizeof(float));
    cudaMalloc((void**)&d_c, n * n * sizeof(float));

    
    cudaMemcpy(d_a, a, n * n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, n * n * sizeof(float), cudaMemcpyHostToDevice);

    const float alpha = 1.0f;
    const float beta = 0.0f;
    cublasHandle_t handle;
    cublasCreate(&handle);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    
    cudaEventRecord(start, 0);
    // Your implementation
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "GPU Elapsed time: " << milliseconds << " ms\n";


    cudaMemcpy(c, d_c, n * n * sizeof(float), cudaMemcpyDeviceToHost);

    // CPU version
    auto start_cpu = std::chrono::high_resolution_clock::now();
    matrixMultiplyCPU(a, b, c_serial, n);
    auto stop_cpu = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> cpu_ms = stop_cpu - start_cpu;
    std::cout << "CPU Elapsed time: " << cpu_ms.count() << " ms\n";

    if (compareMatrices(c, c_serial, n)) {
        std::cout << "The matrices are approximately equal." << std::endl;
    } else {
        std::cout << "There is a discrepancy between the matrices." << std::endl;
    }

    // Cleanup
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);
    free(a);
    free(b);
    free(c);
    free(c_serial);
    cublasDestroy(handle);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return 0;
}
