#include <iostream>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cmath>
#include <string>
#include <climits>

using namespace std;
using namespace std::chrono;

// Utility function to initialize an array with random values
void initializeArray(int* arr, int n, int max_value = 100) {
    for (int i = 0; i < n; i++) {
        arr[i] = rand() % max_value;
    }
}

void initializeArray_prob10(int* arr, int n, int max_value = 999990) {
    for (int i = 0; i < n; i++) {
        arr[i] = (rand() % (2 * max_value + 1)) - max_value;  // Values between -max_value and +max_value
    }
}

// Function to calculate serial and parallel sum
void arraySummation() {
    int n = 1000000;
    int* arr = new int[n];
    for (int i = 0; i < n; i++) arr[i] = i + 1;

    // Serial Sum
    int sum = 0;
    auto start = high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        sum += arr[i];
    }
    auto end = high_resolution_clock::now();
    cout << "Serial Sum: " << sum << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    // Parallel Sum
    sum = 0;
    start = high_resolution_clock::now();
    
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i++) {
        sum += arr[i];
    }

    end = high_resolution_clock::now();
    cout << "Parallel Sum: " << sum << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    delete[] arr;
}

void findMaximumElement() {
    int n = 1000000;
    int* arr = new int[n];
    for (int i = 0; i < n; i++) arr[i] = rand() % 1000000;

    // Serial Max
    int max_val = arr[0];
    auto start = high_resolution_clock::now();
    for (int i = 1; i < n; i++) {
        if (arr[i] > max_val) {
            max_val = arr[i];
        }
    }
    auto end = high_resolution_clock::now();
    cout << "Serial Max: " << max_val << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    // Parallel Max
    max_val = arr[0];
    start = high_resolution_clock::now();
    #pragma omp parallel for reduction(max:max_val)
    for (int i = 1; i < n; i++) {
        if (arr[i] > max_val) {
            max_val = arr[i];
        }
    }
    end = high_resolution_clock::now();
    cout << "Parallel Max: " << max_val << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    delete[] arr;
}

// Function to perform a prime number check
void primeCheck() {
    int number = 10345689;  // A large prime number
    bool is_prime = true;
    int sqrt_number = (int)sqrt(number);
	
    // Serial Prime Check
    auto start = high_resolution_clock::now();
    for (int i = 2; i <= sqrt_number; i++) {
        if (number % i == 0) {
            is_prime = false;
            break;
        }
    }
    auto end = high_resolution_clock::now();
    cout << "Serial Prime Check: " << (is_prime ? "Prime" : "Not Prime") << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    is_prime = true;
    start = high_resolution_clock::now();
    #pragma omp parallel for shared(is_prime)
    for (int i = 2; i <= sqrt_number; i++) {
        if (number % i == 0) {
            #pragma omp critical
            is_prime = false;
        }
    }
    end = high_resolution_clock::now();
    cout << "Parallel Prime Check: " << (is_prime ? "Prime" : "Not Prime") << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";
}

void vectorDotProduct() {
    int n = 1000000;
    int* vec1 = new int[n];
    int* vec2 = new int[n];

    initializeArray(vec1, n);
    initializeArray(vec2, n);

    // Serial Dot Product
    int dot_product = 0;
    auto start = high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        dot_product += vec1[i] * vec2[i];
    }
    auto end = high_resolution_clock::now();
    cout << "Serial Dot Product: " << dot_product << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    // Parallel Dot Product
    dot_product = 0;
    start = high_resolution_clock::now();
    #pragma omp parallel for reduction(+:dot_product)
    for (int i = 0; i < n; i++) {
        dot_product += vec1[i] * vec2[i];
    }
    end = high_resolution_clock::now();
    cout << "Parallel Dot Product: " << dot_product << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    delete[] vec1;
    delete[] vec2;
}

void arrayDifference() {
    int n = 1000000;
    int* arr = new int[n];
    initializeArray(arr, n);

    // Serial Version
    long long diff_serial = arr[0];
    auto start = high_resolution_clock::now();
    for (int i = 1; i < n; i++) {
        diff_serial -= arr[i];
    }
    auto end = high_resolution_clock::now();
    cout << "Serial Cumulative Difference: " << diff_serial 
         << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    // Parallel Version
    long long diff_parallel = arr[0];  
    start = high_resolution_clock::now();
    #pragma omp parallel for reduction(-:diff_parallel)
    for (int i = 1; i < n; i++) {
        diff_parallel -= arr[i];  // Parallel subtraction
    }
    end = high_resolution_clock::now();

    cout << "Parallel Cumulative Difference: " << diff_parallel 
         << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    delete[] arr;
}


void bubbleSort() {
    int n = 10000;  // Size of the array for sorting
    int* arr = new int[n];

    initializeArray(arr, n, 1000);  // Initialize array with values between 0 and 999

    // Serial Bubble Sort
    auto start = high_resolution_clock::now();
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                swap(arr[j], arr[j + 1]);
            }
        }
    }
    auto end = high_resolution_clock::now();
    cout << "Serial Bubble Sort | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    initializeArray(arr, n, 1000);  // Reinitialize array for parallel version

    // Parallel Bubble Sort using OpenMP
    bool sorted = false;
    start = high_resolution_clock::now();
    
    while (!sorted) {
        sorted = true;

        // Odd phase
        #pragma omp parallel for
        for (int i = 1; i < n - 1; i += 2) {
            if (arr[i] > arr[i + 1]) {
                swap(arr[i], arr[i + 1]);
                sorted = false;
            }
        }

        // Even phase
        #pragma omp parallel for
        for (int i = 0; i < n - 1; i += 2) {
            if (arr[i] > arr[i + 1]) {
                swap(arr[i], arr[i + 1]);
                sorted = false;
            }
        }
    }

    end = high_resolution_clock::now();
    cout << "Parallel Bubble Sort | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    delete[] arr;
}

void arraySubtraction() {
    int n = 1000000;  // Size of the arrays
    int* arr1 = new int[n];
    int* arr2 = new int[n];
    int* arr_serial = new int[n];
    int* arr_parallel = new int[n];

    
    initializeArray(arr1, n);
    initializeArray(arr2, n);

    // ==================== Serial Array Subtraction ====================
    auto start = high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        arr_serial[i] = arr1[i] - arr2[i];
    }
    auto end = high_resolution_clock::now();
    cout << "Serial Array Subtraction | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    // Parallel Version
    start = high_resolution_clock::now();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        arr_parallel[i] = arr1[i] - arr2[i];
    }
    end = high_resolution_clock::now();   
    cout << "Parallel Array Subtraction | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    // Deallocate arrays
    delete[] arr1;
    delete[] arr2;
    delete[] arr_serial;
    delete[] arr_parallel;
}

void monteCarloPiEstimation() {
    long long points_inside_circle_serial = 0;
    long long points_inside_circle_parallel = 0;
    int num_points=100000;
   

    // ========== Serial Monte Carlo Pi Estimation ==========
    auto start = high_resolution_clock::now();
    for (long long i = 0; i < num_points; i++) {
        // Generate random x and y between -1 and 1
        double x = (double)rand() / RAND_MAX * 2.0 - 1.0;
        double y = (double)rand() / RAND_MAX * 2.0 - 1.0;

        // Check if the point is inside the circle
        if (x * x + y * y <= 1) {
            points_inside_circle_serial++;
        }
    }
    auto end = high_resolution_clock::now();
    double pi_serial = 4.0 * points_inside_circle_serial / num_points;
    cout << "Serial Pi Estimation: " << pi_serial << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    // ========== Parallel Monte Carlo Pi Estimation Using OpenMP ==========
    points_inside_circle_parallel = 0;  // Reset counter
    start = high_resolution_clock::now();

    #pragma omp parallel
    {
        // Thread-local counter
        long long local_count = 0;
        unsigned int seed = time(0) ^ omp_get_thread_num();  // Unique seed for each thread

        #pragma omp for
        for (long long i = 0; i < num_points; i++) {
            // Generate random x and y using thread-local seed
            double x = (double)rand_r(&seed) / RAND_MAX * 2.0 - 1.0;
            double y = (double)rand_r(&seed) / RAND_MAX * 2.0 - 1.0;

            if (x * x + y * y <= 1) {
                local_count++;
            }
        }

        // Accumulate local counts to the global count
        #pragma omp atomic
        points_inside_circle_parallel += local_count;
    }
   
    end = high_resolution_clock::now();
    double pi_parallel = 4.0 * points_inside_circle_parallel / num_points;
    cout << "Parallel Pi Estimation: " << pi_parallel << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";
}

void sumOfSquaresCustomReduction() {
    int n=1000000;
    int* arr_s = new int[n];
    initializeArray(arr_s, n);

    // ==================== Serial Sum of Squares ====================
    int sum_serial = 0;
    auto start = high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        sum_serial += arr_s[i] * arr_s[i];
    }
    auto end = high_resolution_clock::now();
    cout << "Serial Sum of Squares: " << sum_serial << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    // ==================== Parallel Sum of Squares without build in OpenMP reduction clause (You can use any other directives)) ====================
   
    start = high_resolution_clock::now();
    
    int sum_builtin = 0;
    start = high_resolution_clock::now();
    #pragma omp parallel for reduction(+:sum_builtin)
    for (int i = 0; i < n; i++) {
        sum_builtin += arr_s[i] * arr_s[i];
    }
    end = high_resolution_clock::now();
    cout << "Parallel Sum of Squares (Custom Reduction): " << sum_builtin << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

}

void findMinimumValue() {
    int n = 1000000;
    int* arr = new int[n];
    initializeArray_prob10(arr, n);
    
    // Serial Code for Finding the Minimum Value
    int min_value_serial = INT_MAX;
    auto start = high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
        if (arr[i] < min_value_serial) {
            min_value_serial = arr[i];
        }
    }
    auto end = high_resolution_clock::now();
    cout << "Serial Minimum Value: " << min_value_serial << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";
    
    
    // Parallel Code for Finding the Minimum Value using Custom Reduction function clause
    // Note: If you want full points for this problem, you will need implement the solution with Custom Reduction function
    // Hint: Use #pragma omp declare reduction.
    int min_value_parallel = INT_MAX;
    start = high_resolution_clock::now();

    #pragma omp declare reduction(minimum:int:omp_out = omp_in < omp_out ? omp_in : omp_out) initializer(omp_priv = INT_MAX)
    
    // Parallel Code for Finding the Minimum Value using Custom Reduction
    #pragma omp parallel for reduction(minimum:min_value_parallel)
    for (int i = 0; i < n; i++) {
        if (arr[i] < min_value_parallel) {
            min_value_parallel = arr[i];
        }
    }

    end = high_resolution_clock::now();
    cout << "Parallel Minimum Value: " << min_value_parallel << " | Time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";
    
   
    delete[] arr;
}

int main() {
    srand(time(0));  // Seed for random number generation

    cout << "\n================= Problem 1: Array Summation =================\n";
    arraySummation();

    cout << "\n================= Problem 2: Finding the Maximum Element =================\n";
    findMaximumElement();

    cout << "\n================= Problem 3: Prime Number Check =================\n";
    primeCheck();

    cout << "\n================= Problem 4: Vector Dot Product =================\n";
    vectorDotProduct();

    cout << "\n================= Problem 5: Difference of elements in array =================\n";
    arrayDifference();

    cout << "\n================= Problem 6: Bubble Sort =================\n";
    bubbleSort();

    cout << "\n================= Problem 7: Array Substraction =================\n";
    arraySubtraction();

    cout << "\n================= Problem 8: Monte Carlo Estimation  =================\n";
    monteCarloPiEstimation();

    cout << "\n================= Problem 9: Sum of Squares  =================\n";
    sumOfSquaresCustomReduction();

    cout << "\n================= Problem 10: Finding Minimum Value (Using Reduction Function)  =================\n";
    findMinimumValue();

    return 0;
}