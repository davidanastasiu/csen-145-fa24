#include <cmath>
#include <vector>
#include <numeric>
#include <random>
#include <algorithm>
#include <execution>
#include <iostream>

int main(void) {
    size_t n = 10000;

    // Generate a vector of length n
    std::vector<double> v(n);

    // Initilaize the vector with random values
    // Initial random values were too large and caused overflow
    std::mt19937 engine{std::random_device{}()};
    std::uniform_int_distribution dist{1, 100}; // inclusive
    std::generate(v.begin(), v.end(), [&] { return dist(engine); });

    // Compute the sum of the numbers
    int sum = std::reduce(std::execution::par, v.begin(), v.end(), 0);
    std::cout << "Sum of v: " << sum << std::endl;
    
    // Make a copy of v
    auto v2 = v;

    // square the vector v via foreach
    // Create a vector of indicies
    std::vector<size_t> i = std::vector<size_t>(n);
    // std::iota(std::begin(i), std::end(i), 0);
    std::generate(std::execution::par, i.begin(), i.end(), [n = 0] () mutable { return n++; });

    // print the first 10 elements of the vector
    for(int i=0; i < 10; ++i){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;

    // Square the vector v via for_each
    std::for_each(std::execution::par, i.begin(), i.end(), [&](auto& j) {
        v[j] = v[j] * v[j]; 
    });

    // print the first 10 elements of the vector
    for(int i=0; i < 10; ++i){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;

    // std::execution::seq: Sequential execution.
    // std::execution::par: Parallel execution.
    // std::execution::par_unseq: Parallel and potentially vectorized execution.

    // print the sum of the squared vector v
    sum = std::reduce(std::execution::par, v.begin(), v.end(), 0);
    std::cout << "Sum of v^2: " << sum << std::endl;

    // square the vector v2 via transform
    std::transform(std::execution::par, v2.begin(), v2.end(), v2.begin(),
    [](int num) { 
        return num * num; 
    });

    // print the first 10 elements of the vector
    for(int i=0; i < 10; ++i){
        std::cout << v2[i] << " ";
    }
    std::cout << std::endl;

    // print the sum of the squared vector v2
    sum = std::reduce(std::execution::par, v2.begin(), v2.end(), 0);
    std::cout << "Sum of v2^2: " << sum << std::endl;

    return EXIT_SUCCESS;
}