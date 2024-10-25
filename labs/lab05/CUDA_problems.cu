#include <cuda_runtime.h>
#include <iostream>
#include <chrono>
#include <cstdlib>

#define N 1000000  // Array size for vector operations
#define BLOCK_SIZE 256     // Array size for vector operations
#define WIDTH 1024     // Width of matrix/image
#define HEIGHT 1024    // Height of matrix/image 

// Utility function to check for CUDA errors
#define CUDA_CHECK(call) { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error in " << __FILE__ << " at line " << __LINE__ << ": " \
                  << cudaGetErrorString(err) << std::endl; \
        exit(EXIT_FAILURE); \
    } \
}

// 1. Vector Addition Kernel
__global__ void vectorAdd(float* a, float* b, float* c, int n) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    // Your Implementation
}

void vectorAddCPU(float* a, float* b, float* c, int n) {
    for (int i = 0; i < n; ++i) {
        c[i] = a[i] + b[i];
    }
}


__global__ void histogramKernel(int* input, int* histogram, int n, int bins) {
    // Your Implementation
    // Hint: Use AtomicAdd
}

void histogramCPU(int* input, int* histogram, int n, int bins) {
    for (int i = 0; i < bins; ++i) histogram[i] = 0;
    for (int i = 0; i < n; ++i) histogram[input[i]]++;
}

__global__ void nearestNeighborKernel(float* data, float* query, float* result, int n) {
    // Your Implementation
     // Hint: Use AtomicMin
}

void nearestNeighborCPU(float* data, float query, float& result, int n) {
    result = abs(data[0] - query);
    for (int i = 1; i < n; ++i) {
        float diff = abs(data[i] - query);
        if (diff < result) result = diff;
    }
}



__global__ void reverseArray(float* input, float* output, int n) {
    // Your Implementation
}

void reverseArrayCPU(float* input, float* output, int n) {
    for (int i = 0; i < n; ++i) {
        output[n - i - 1] = input[i];
    }
}

__global__ void transposeKernel(float* input, float* output, int width, int height) {
    
    // Your Implementation
}

void transposeCPU(float* input, float* output, int width, int height) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            output[j * height + i] = input[i * width + j];
        }
    }
}

__global__ void convolutionKernel(float* image, float* output, float* filter, int width, int height) {
    // Your Implementation
}

// Convolution (CPU)
void convolutionCPU(float* image, float* output, float* filter, int width, int height) {
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float sum = 0.0f;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    int nx = std::min(std::max(x + i, 0), width - 1);
                    int ny = std::min(std::max(y + j, 0), height - 1);
                    sum += image[ny * width + nx] * filter[(i + 1) * 3 + (j + 1)];
                }
            }
            output[y * width + x] = sum;
        }
    }
}


// 7. Bitwise AND Operation
__global__ void bitwiseAnd(int* a, int* b, int* c, int n) {
    // Your Implementation
}

void bitwiseAndCPU(int* a, int* b, int* c, int n) {
    for (int i = 0; i < n; ++i) {
        c[i] = a[i] & b[i];
    }
}
// Utility function to measure time
template <typename Func>
void measureTime(Func func, const char* description) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> duration = end - start;
    std::cout << description << " took " << duration.count() << " ms\n";
}

template <typename T>
void validateResults(T* gpuResult, T* cpuResult, int n, const std::string& operationName) {
    for (int i = 0; i < n; ++i) {
        if (std::abs(gpuResult[i] - cpuResult[i]) > 1e-5) {
            std::cerr << operationName << " failed at index " << i 
                      << ": GPU result = " << gpuResult[i] 
                      << ", CPU result = " << cpuResult[i] << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    std::cout << operationName << " results match!\n";
}

int main() {
    // Allocate and initialize host memory
    float *h_a = new float[N], *h_b = new float[N], *h_cGPU = new float[N], *h_cCPU = new float[N];
    int *h_int_a = new int[N], *h_int_b = new int[N], *h_int_cGPU = new int[N], *h_int_cCPU = new int[N];
    
    for (int i = 0; i < N; ++i) {
        h_a[i] = h_b[i] = static_cast<float>(rand()) / RAND_MAX;
        h_int_a[i] = h_int_b[i] = rand() % 256;
    }

    // Allocate device memory
    int *h_input = new int[N], *h_histogramCPU = new int[256], *h_histogramGPU = new int[256];
    int *d_input, *d_histogram;
    float *d_a, *d_b, *d_c;
    int *d_int_a, *d_int_b, *d_int_c;
    
    CUDA_CHECK(cudaMalloc(&d_a, N * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_b, N * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_c, N * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_int_a, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_int_b, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_int_c, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_input, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_histogram, 256 * sizeof(int)));
    // -------------------------------------------------------------------------------------------------------------------------//
    // Initialize input for histogram
    for (int i = 0; i < N; ++i) h_input[i] = rand() % 256;
    CUDA_CHECK(cudaMemcpy(d_input, h_input, N * sizeof(int), cudaMemcpyHostToDevice));

    // GPU Histogram
    measureTime([&]() {
        CUDA_CHECK(cudaMemset(d_histogram, 0, 256 * sizeof(int)));
        histogramKernel<<<(N + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>(d_input, d_histogram, N, 256);
        CUDA_CHECK(cudaDeviceSynchronize());
    }, "GPU Histogram");

    CUDA_CHECK(cudaMemcpy(h_histogramGPU, d_histogram, 256 * sizeof(int), cudaMemcpyDeviceToHost));

    // CPU Histogram
    measureTime([&]() { histogramCPU(h_input, h_histogramCPU, N, 256); }, "CPU Histogram");
    
    // Validate Histogram Results
    validateResults(h_histogramGPU, h_histogramCPU, 256, "Histogram");

    // Copy data from host to device for vector operations
    CUDA_CHECK(cudaMemcpy(d_a, h_a, N * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, N * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_int_a, h_int_a, N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_int_b, h_int_b, N * sizeof(int), cudaMemcpyHostToDevice));

    // -------------------------------------------------------------------------------------------------------------------------//

    // Vector Addition
    measureTime([&]() {
        vectorAdd<<<(N + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>(d_a, d_b, d_c, N);
        CUDA_CHECK(cudaDeviceSynchronize());
    }, "GPU Vector Addition");

    CUDA_CHECK(cudaMemcpy(h_cGPU, d_c, N * sizeof(float), cudaMemcpyDeviceToHost));

    measureTime([&]() { vectorAddCPU(h_a, h_b, h_cCPU, N); }, "CPU Vector Addition");
    validateResults(h_cGPU, h_cCPU, N, "Vector Addition");

    // -------------------------------------------------------------------------------------------------------------------------//
    // Bitwise AND
    measureTime([&]() {
        bitwiseAnd<<<(N + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>(d_int_a, d_int_b, d_int_c, N);
        CUDA_CHECK(cudaDeviceSynchronize());
    }, "GPU Bitwise AND");

    CUDA_CHECK(cudaMemcpy(h_int_cGPU, d_int_c, N * sizeof(int), cudaMemcpyDeviceToHost));

    measureTime([&]() { bitwiseAndCPU(h_int_a, h_int_b, h_int_cCPU, N); }, "CPU Bitwise AND");
    validateResults(h_int_cGPU, h_int_cCPU, N, "Bitwise AND");


    float *h_reverseInput = new float[N];  // Host input array for reverse
    float *h_reverseOutputGPU = new float[N];  // Host GPU result
    float *h_reverseOutputCPU = new float[N];  // Host CPU result

    // Initialize input array with random values
    for (int i = 0; i < N; ++i) {
        h_reverseInput[i] = static_cast<float>(rand()) / RAND_MAX;
    }

    // Allocate device memory for Array Reverse
    float *d_reverseInput, *d_reverseOutput;
    CUDA_CHECK(cudaMalloc(&d_reverseInput, N * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_reverseOutput, N * sizeof(float)));

    // Copy input data to the device
    CUDA_CHECK(cudaMemcpy(d_reverseInput, h_reverseInput, N * sizeof(float), cudaMemcpyHostToDevice));

    // Define grid and block dimensions
    // Execute the Array Reverse on GPU
    measureTime([&]() {
        reverseArray<<<(N + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE>>>(d_reverseInput, d_reverseOutput, N);
        CUDA_CHECK(cudaDeviceSynchronize());
    }, "GPU Array Reverse");

    // Copy the result back to the host
    CUDA_CHECK(cudaMemcpy(h_reverseOutputGPU, d_reverseOutput, N * sizeof(float), cudaMemcpyDeviceToHost));

    // Execute the Array Reverse on CPU
    measureTime([&]() { reverseArrayCPU(h_reverseInput, h_reverseOutputCPU, N); }, "CPU Array Reverse");

    // Validate the Array Reverse Results
    validateResults(h_reverseOutputGPU, h_reverseOutputCPU, N, "Array Reverse");

   // -------------------------------------------------------------------------------------------------------------------------//

    // Nearest Neighbor search
    float query = static_cast<float>(rand()) / RAND_MAX;  // Random query value
    float nearestNeighborResultGPU = std::numeric_limits<float>::max();
    float nearestNeighborResultCPU = std::numeric_limits<float>::max();
    float *d_data, *d_query, *d_result;
    CUDA_CHECK(cudaMalloc(&d_data, N * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_query, sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_result, sizeof(float)));

    // Copy data and query to the device
    CUDA_CHECK(cudaMemcpy(d_data, h_a, N * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_query, &query, sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_result, &nearestNeighborResultGPU, sizeof(float), cudaMemcpyHostToDevice));


    // Execute the Nearest Neighbor Search on GPU
    dim3 blockSize(BLOCK_SIZE);
    dim3 gridSize((N + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // Execute the Nearest Neighbor Search on GPU
    measureTime([&]() {
        nearestNeighborKernel<<<gridSize, blockSize>>>(d_data, d_query, d_result, N);
        CUDA_CHECK(cudaDeviceSynchronize());
    }, "GPU Nearest Neighbor Search");
    // Copy the result back to the host
    CUDA_CHECK(cudaMemcpy(&nearestNeighborResultGPU, d_result, sizeof(float), cudaMemcpyDeviceToHost));

    // Execute the Nearest Neighbor Search on CPU
    measureTime([&]() { nearestNeighborCPU(h_a, query, nearestNeighborResultCPU, N); }, 
                "CPU Nearest Neighbor Search");

    // Validate the Nearest Neighbor Search Results
    if (std::abs(nearestNeighborResultGPU - nearestNeighborResultCPU) > 1e-5) {
        std::cerr << "Nearest Neighbor Search failed: GPU result = " 
                  << nearestNeighborResultGPU << ", CPU result = " 
                  << nearestNeighborResultCPU << std::endl;
    } else {
        std::cout << "Nearest Neighbor Search results match!\n";
    }

    float *h_image = new float[WIDTH * HEIGHT], *h_outputGPU = new float[WIDTH * HEIGHT], *h_outputCPU = new float[WIDTH * HEIGHT];
    float filter[9] = {0, -1, 0, -1, 5, -1, 0, -1, 0};  // Simple sharpening filter

    float *d_image, *d_output, *d_filter;
    CUDA_CHECK(cudaMalloc(&d_image, WIDTH * HEIGHT * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_output, WIDTH * HEIGHT * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_filter, 9 * sizeof(float)));

    // Initialize the image with random values
    for (int i = 0; i < WIDTH * HEIGHT; ++i) h_image[i] = static_cast<float>(rand()) / RAND_MAX;

    // Copy image and filter to device memory
    CUDA_CHECK(cudaMemcpy(d_image, h_image, WIDTH * HEIGHT * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_filter, filter, 9 * sizeof(float), cudaMemcpyHostToDevice));

    // === Convolution Execution and Validation ===
    measureTime([&]() {
        dim3 blockSize(16, 16);
        dim3 gridSize((WIDTH + blockSize.x - 1) / blockSize.x, (HEIGHT + blockSize.y - 1) / blockSize.y);
        convolutionKernel<<<gridSize, blockSize>>>(d_image, d_output, d_filter, WIDTH, HEIGHT);
        CUDA_CHECK(cudaDeviceSynchronize());
    }, "GPU Convolution");

    // Copy result back to host
    CUDA_CHECK(cudaMemcpy(h_outputGPU, d_output, WIDTH * HEIGHT * sizeof(float), cudaMemcpyDeviceToHost));

    // Run convolution on CPU
    measureTime([&]() { convolutionCPU(h_image, h_outputCPU, filter, WIDTH, HEIGHT); }, "CPU Convolution");

    // Validate results
    validateResults(h_outputGPU, h_outputCPU, WIDTH * HEIGHT, "Convolution");

    // -------------------------------------------------------------------------------------------------------------------------//
    // === Matrix Transposition Setup ===
    float *h_matrix = new float[WIDTH * HEIGHT], *h_transposedGPU = new float[WIDTH * HEIGHT], *h_transposedCPU = new float[WIDTH * HEIGHT];

    float *d_matrix, *d_transposed;
    CUDA_CHECK(cudaMalloc(&d_matrix, WIDTH * HEIGHT * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_transposed, WIDTH * HEIGHT * sizeof(float)));

    // Initialize the matrix with random values
    for (int i = 0; i < WIDTH * HEIGHT; ++i) h_matrix[i] = static_cast<float>(rand()) / RAND_MAX;

    // Copy matrix to device memory
    CUDA_CHECK(cudaMemcpy(d_matrix, h_matrix, WIDTH * HEIGHT * sizeof(float), cudaMemcpyHostToDevice));

    // === Matrix Transposition Execution and Validation ===
    measureTime([&]() {
        dim3 blockSize(16, 16);
        dim3 gridSize((WIDTH + blockSize.x - 1) / blockSize.x, (HEIGHT + blockSize.y - 1) / blockSize.y);
        transposeKernel<<<gridSize, blockSize>>>(d_matrix, d_transposed, WIDTH, HEIGHT);
        CUDA_CHECK(cudaDeviceSynchronize());
    }, "GPU Matrix Transposition");

    // Copy result back to host
    CUDA_CHECK(cudaMemcpy(h_transposedGPU, d_transposed, WIDTH * HEIGHT * sizeof(float), cudaMemcpyDeviceToHost));

    // Run transposition on CPU
    measureTime([&]() { transposeCPU(h_matrix, h_transposedCPU, WIDTH, HEIGHT); }, "CPU Matrix Transposition");

    // Validate results
    validateResults(h_transposedGPU, h_transposedCPU, WIDTH * HEIGHT, "Matrix Transposition");
    

    // Cleanup
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_c));
    CUDA_CHECK(cudaFree(d_int_a));
    CUDA_CHECK(cudaFree(d_int_b));
    CUDA_CHECK(cudaFree(d_int_c));
    CUDA_CHECK(cudaFree(d_input));
    CUDA_CHECK(cudaFree(d_histogram));
    CUDA_CHECK(cudaFree(d_image));
    CUDA_CHECK(cudaFree(d_transposed));
    CUDA_CHECK(cudaFree(d_reverseInput));
    CUDA_CHECK(cudaFree(d_reverseOutput));
    CUDA_CHECK(cudaFree(d_output));
    CUDA_CHECK(cudaFree(d_filter));
    CUDA_CHECK(cudaFree(d_matrix));
    CUDA_CHECK(cudaFree(d_data));
    CUDA_CHECK(cudaFree(d_query));
    CUDA_CHECK(cudaFree(d_result));

    delete[] h_a;
    delete[] h_b;
    delete[] h_cGPU;
    delete[] h_cCPU;
    delete[] h_int_a;
    delete[] h_int_b;
    delete[] h_int_cGPU;
    delete[] h_int_cCPU;
    delete[] h_input;
    delete[] h_histogramCPU;
    delete[] h_histogramGPU;
    delete[] h_image;
    delete[] h_outputGPU;
    delete[] h_outputCPU;
    delete[] h_matrix;
    delete[] h_transposedGPU;
    delete[] h_transposedCPU;
    delete[] h_reverseInput;
    delete[] h_reverseOutputGPU;
    delete[] h_reverseOutputCPU;

    return 0;
}
