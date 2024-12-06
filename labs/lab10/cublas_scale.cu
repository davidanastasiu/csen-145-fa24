#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#define M 1000 
#define N 700

// Function to scale matrix serially
void scaleMatrixSerial(std::vector<float>& matrix, float alpha, int rows, int cols) {
    for (int i = 0; i < rows * cols; ++i) {
        matrix[i] *= alpha;
    }
}


bool areMatricesEqual(const std::vector<float>& a, const std::vector<float>& b, float tolerance) {
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i] - b[i]) > tolerance) {
            return false;
        }
    }
    return true;
}

int main() {
    std::vector<float> a(M * N), b(M * N);
    float* devPtrA = nullptr;
    cublasHandle_t handle;

    // Factor by which to scale the matrix
    float alpha = 10.0f;

   
    srand(42); 

    for (int i = 0; i < M * N; ++i) {
        float randomValue = static_cast<float>(rand()) / RAND_MAX; // Generate a random float between 0 and 1
        a[i] = randomValue;
        b[i] = randomValue; 
    }
        cudaEvent_t startCUDA, stopCUDA;
        cudaEventCreate(&startCUDA);
        cudaEventCreate(&stopCUDA);
        float millisecondsCUDA = 0;
    try {
        cudaMalloc((void**)&devPtrA, M * N * sizeof(float));
        cublasCreate(&handle);
        cublasSetMatrix(M, N, sizeof(float), a.data(), M, devPtrA, M);
        cudaEventRecord(startCUDA);
        // Your implementation
        cudaEventRecord(stopCUDA);
        cudaEventSynchronize(stopCUDA);
        cudaEventElapsedTime(&millisecondsCUDA, startCUDA, stopCUDA);
        std::cout << "GPU Processing time: " << millisecondsCUDA << " ms\n";
        cublasGetMatrix(M, N, sizeof(float), devPtrA, M, a.data(), M);
        auto startSerial = std::chrono::high_resolution_clock::now();
        scaleMatrixSerial(b, alpha, M, N);
        auto endSerial = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> millisecondsSerial = endSerial - startSerial;
        std::cout << "Serial computation time: " << millisecondsSerial.count() << " ms\n";

        if (areMatricesEqual(a, b, 1e-5)) {
            std::cout << "The matrices are approximately equal." << std::endl;
        } else {
            std::cout << "There is a discrepancy between the matrices." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "An exception occurred: " << e.what() << std::endl;
        if (devPtrA) cudaFree(devPtrA);
        if (handle) cublasDestroy(handle);
        return EXIT_FAILURE;
    }

    
    cudaFree(devPtrA);
    cublasDestroy(handle);
    cudaEventDestroy(startCUDA);
    cudaEventDestroy(stopCUDA);

    return EXIT_SUCCESS;
}
