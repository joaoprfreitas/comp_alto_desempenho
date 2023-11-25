all:
	mpicc main.c -o par -fopenmp -lm -Wall
	gcc dist-seq.c -o seq -lm

mpi:
	mpirun -np 4 ./par 30 1 6

test:
	./seq

clean:
	rm seq par