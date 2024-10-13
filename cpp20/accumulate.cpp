#include <iostream>
#include <vector>
#include <numeric>
#include <execution>

int main(void) {
  size_t n = 10000;

  // Generate a vector of length n
  std::vector<int> v = std::vector<int>(n);

  // Initilaize the vector with values
  std::fill(std::execution::par, v.begin(), v.end(), -1);

  // Compute the sum of all elements
  int s = std::reduce(std::execution::par, v.begin(), v.end(), 0.0);
  
  // Output the result
  std::cout << "Result= " << s << std::endl;

  return EXIT_SUCCESS;
}