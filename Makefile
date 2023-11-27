all:
	mpicc projeto.c -o par -fopenmp -lm -Wall
	gcc dist-seq.c -o seq -lm

mpi:
	mpicc projeto.c -o par -fopenmp -lm -Wall
	mpirun -np 4 ./par 10 1 4

test:
	gcc dist-seq.c -o seq -lm
	./seq 100 1

compare:
	gcc dist-seq.c -o seq -lm
	mpicc main6.c -o par -fopenmp -lm -Wall
	./seq
	mpirun -np 4 ./par 100 1 2

clean:
	rm seq par