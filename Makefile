mpi:
	mpicc main.c -o main -fopenmp -lm -Wall
	mpirun --oversubscribe -np 2 ./main 100 1 4

test:
	./testecase