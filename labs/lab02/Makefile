CC = g++

default: DenseMatMult

DenseMatMult: DenseMatMult.cpp
	${CC} -O3 -Wall -Wextra -std=c++20 -o $@ DenseMatMult.cpp

DenseMatMult2: DenseMatMult2.cpp
	${CC} -O3 -Wall -Wextra -std=c++20 -o $@ DenseMatMult.cpp

DenseMatMultOmp: DenseMatMultOmp.cpp
	${CC} -O3 -Wall -Wextra -std=c++20 -fopenmp -o $@ DenseMatMultOmp.cpp

DenseMatMult2Omp: DenseMatMult2Omp.cpp
	${CC} -O3 -Wall -Wextra -std=c++20 -fopenmp -o $@ DenseMatMult2Omp.cpp

clean:
	-rm -f DenseMatMult DenseMatMult2 DenseMatMultOmp DenseMatMult2Omp

