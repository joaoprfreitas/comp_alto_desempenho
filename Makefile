all:
	mpicc main.c -o main -fopenmp -lm -Wall

mpi:
	mpirun --oversubscribe -np 4 ./main 200 1 10

test:
	gcc dist-seq-teste.c -o teste -lm
	./teste

clean:
	rm main teste

time:
	mpicc main.c -o main -fopenmp -lm -Wall
	gcc dist-seq-teste.c -o teste -lm
	./teste
	mpirun --oversubscribe -np 4 ./main 100 1 5