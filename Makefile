mpi:
	mpicc main.c -o main -fopenmp -lm -Wall
	mpirun --oversubscribe -np 2 ./main 100 1 10

valgrind:
	mpicc main.c -o main -fopenmp -lm -Wall
	mpirun --oversubscribe -np 1 valgrind --tool=memcheck --leak-check=full ./main 100 1 10

test:
	gcc dist-seq-test.c -o teste.c -lm
	./testecase