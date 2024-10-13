CC = g++

default: accumulate foreach

accumulate: accumulate.cpp
	${CC} -O3 -Wall -Wextra -std=c++20 -o $@ accumulate.cpp

foreach: foreach.cpp
	${CC} -O3 -Wall -Wextra -std=c++20 -o $@ foreach.cpp

clean:
	-rm -f accumulate foreach

