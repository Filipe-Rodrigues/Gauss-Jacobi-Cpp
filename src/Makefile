# Makefile - jacobi #

all: gauss_jacobi

gauss_jacobi: main.cpp gauss_jacobi.cpp
	g++ -std=gnu++11 -fopenmp -pthread main.cpp -o jacobi

clean:
	rm *.o gauss_jacobi