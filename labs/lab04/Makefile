CC = g++

default: sparsematmult

sparsematmult: sparsematmult.cpp
	${CC} -g -Wall -Wextra -Wno-unused-parameter -fopenmp -o $@ sparsematmult.cpp

clean:
	-rm -f sparsematmult
