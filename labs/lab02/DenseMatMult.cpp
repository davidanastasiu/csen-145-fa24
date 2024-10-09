#include <iostream>
#include <chrono>  // For measuring execution time
#include <cstdlib> // For rand() and srand()
#include <ctime>   
using namespace std;
using namespace std::chrono;

// Function to allocate memory for a matrix
int** allocateMatrix(int rows, int cols) {
    int** matrix = new int*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new int[cols];
    }
    return matrix;
}

// Function to fill matrix with random values between 1 and 10
void fillMatrixWithRandomValues(int** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = rand() % 10 + 1;  
        }
    }
}

// Function to multiply two matrices
int** multiplyMatrices(int** mat1, int** mat2, int rows1, int cols1, int cols2) {
    
    auto mat3 = allocateMatrix(rows1, cols2);
    for(int i=0; i < rows1; ++i){
        for(int j=0; j < cols2; ++j){
            int v = 0;
            for(int k=0; k < cols1; ++k){
                v += mat1[i][k] * mat2[k][j];
            }
            mat3[i][j] = v;
        }
    }
    return mat3;
}

// Function to deallocate memory for a matrix
void deallocateMatrix(int** matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

int main() {
    srand(time(0));

    // Test cases with diverse dimensions
    int testCases[10][4] = {
        {348, 556, 556, 648},
        {600, 800, 800, 700},
        {450, 600, 600, 750},
        {512, 512, 512, 600},
        {700, 800, 800, 900},
        {333, 444, 444, 555},
        {688, 432, 432, 532},
        {760, 540, 540, 650},
        {640, 480, 480, 720},
        {950, 850, 850, 950}
    };

    for (int i = 0; i < 10; i++) {
        int rows1 = testCases[i][0];
        int cols1 = testCases[i][1];
        int rows2 = testCases[i][2];
        int cols2 = testCases[i][3];

        // Allocate matrices
        int** mat1 = allocateMatrix(rows1, cols1);
        int** mat2 = allocateMatrix(rows2, cols2);

        // Fill matrices with random values
        fillMatrixWithRandomValues(mat1, rows1, cols1);
        fillMatrixWithRandomValues(mat2, rows2, cols2);

        // Measure the time taken for matrix multiplication
        auto start = high_resolution_clock::now();
        int** result = multiplyMatrices(mat1, mat2, rows1, cols1, cols2);
        auto stop = high_resolution_clock::now();

        // Calculate duration in microseconds
        auto duration = duration_cast<microseconds>(stop - start);

        
        cout << "Test Case " << i + 1 << " (Dimensions: " << rows1 << "x" << cols1 << " * " << rows2 << "x" << cols2 << ") - Time taken for matrix multiplication: " << duration.count() << " microseconds\n";

        // Deallocate memory
        deallocateMatrix(mat1, rows1);
        deallocateMatrix(mat2, rows2);
        deallocateMatrix(result, rows1);
    }

    return 0;
}
