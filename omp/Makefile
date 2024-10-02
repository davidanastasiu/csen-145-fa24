CC = g++

default: examples

examples: examples.cpp
	${CC} -O3 -Wall -Wextra -std=c++20 -fopenmp -o $@ examples.cpp

clean:
	-rm -f examples

