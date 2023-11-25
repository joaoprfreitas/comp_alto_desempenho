mpi:
	mpicc main.c -o main -fopenmp -lm -Wall
	mpirun --oversubscribe -np 4 ./main 10 1 10

test:
	gcc dist-seq-teste.c -o teste -lm
	./teste

time:
	mpicc main.c -o main -fopenmp -lm -Wall
	gcc dist-seq-teste.c -o teste -lm
	./teste
	mpirun --oversubscribe -np 4 ./main 100 1 5