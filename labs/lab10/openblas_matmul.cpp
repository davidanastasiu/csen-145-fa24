#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <cblas.h>

const int SIZE = 1000;


void naiveMultiplication(const std::vector<float>& A, const std::vector<float>& B, std::vector<float>& C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
}
bool compareMatrices(const std::vector<float>& A, const std::vector<float>& B, float tolerance, int N) {
    for (int i = 0; i < N * N; i++) {
        if (fabs(A[i] - B[i]) > tolerance) {
            return false;
        }
    }
    return true;
}

int main() {
    std::vector<float> A(SIZE * SIZE);
    std::vector<float> B(SIZE * SIZE);
    std::vector<float> C(SIZE * SIZE, 0);
    std::vector<float> D(SIZE * SIZE, 0); 

    
    srand(42);
    for (int i = 0; i < SIZE * SIZE; i++) {
        A[i] = static_cast<float>(rand()) / RAND_MAX;
        B[i] = static_cast<float>(rand()) / RAND_MAX;
    }


    auto start = std::chrono::high_resolution_clock::now();
    // Your Implementation
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "OpenBLAS multiplication took: " << elapsed.count() << " seconds." << std::endl;

    
    start = std::chrono::high_resolution_clock::now();
    naiveMultiplication(A, B, D, SIZE);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Naive multiplication took: " << elapsed.count() << " seconds." << std::endl;

    if (compareMatrices(C, D, 1e-3, SIZE)) {
        std::cout << "Results are correct!" << std::endl;
    } else {
        std::cout << "Results differ!" << std::endl;
    }


    return 0;
}
