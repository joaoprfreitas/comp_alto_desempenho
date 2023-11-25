all:
	mpicc main.c -o main -fopenmp -lm -Wall
	gcc dist-seq-teste.c -o teste -lm

mpi:
	mpirun -np 4 ./main 100 1 5

test:
	./teste

clean:
	rm main teste