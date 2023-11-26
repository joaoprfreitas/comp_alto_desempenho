all:
	mpicc main2.c -o par -fopenmp -lm -Wall
	gcc dist-seq.c -o seq -lm

mpi:
	mpirun -np 4 ./par 17 1 2

test:
	./seq

clean:
	rm seq par