#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cblas.h>

double dotProductSerial(const std::vector<double>& x, const std::vector<double>& y) {
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        sum += x[i] * y[i];
    }
    return sum;
}

int main() {
    const int N = 1000000; 
    std::vector<double> x(N), y(N);

    srand(42);

   
    for (int i = 0; i < N; i++) {
        x[i] = static_cast<double>(rand()) / RAND_MAX;
        y[i] = static_cast<double>(rand()) / RAND_MAX;
    }


    auto start = std::chrono::high_resolution_clock::now();
    // Your implementation
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedOpenBLAS = end - start;
    std::cout << "Dot product using OpenBLAS: " << resultOpenBLAS << std::endl;
    std::cout << "OpenBLAS computation time: " << elapsedOpenBLAS.count() << " seconds" << std::endl;

    
    start = std::chrono::high_resolution_clock::now();
    double resultSerial = dotProductSerial(x, y);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedSerial = end - start;
    std::cout << "Dot product using serial implementation: " << resultSerial << std::endl;
    std::cout << "Serial computation time: " << elapsedSerial.count() << " seconds" << std::endl;


    if (std::abs(resultOpenBLAS - resultSerial) < 1e-5) {
        std::cout << "Results are consistent in both implementations." << std::endl;
    } else {
        std::cout << "Results differ between implementations!" << std::endl;
    }

    return 0;
}
