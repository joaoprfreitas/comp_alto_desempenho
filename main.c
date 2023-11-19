/* Integrantes:
 * João Pedro Rodrigues Freitas - 11316552
 * Gabriel Akio Urakawa - 11795912
 * Renato Tadeu Theodoro Junior - 11796750
 * Samuel Victorio Bernacci - 12703455
 * 
 * Comandos para compilar:
 * 
 * 
 * 
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>
#include <limits.h>

#define MAX_VALUE 100

typedef unsigned long long int ulli;
typedef long double ld;

// Retorna a distância euclidiana entre dois pontos (x, y, z)
double euclidian(int x1, int y1, int z1, int x2, int y2, int z2) {
    int x = x1 - x2;
    int y = y1 - y2;
    int z = z1 - z2;
    return sqrt(x*x + y*y + z*z);
}

// Retorna a distância de manhattan entre dois pontos (x, y, z)
int manhattan(int x1, int y1, int z1, int x2, int y2, int z2) {
    return abs(x1 - x2) + abs(y1 - y2) + abs(z1 - z2);
}

int *createMatrix(int row, int column, int SEED);
void printMatrix(int *matrix, int row, int column);
void freeMatrix(int *matrix, int row);

// TODO: tratar o caso de N n caber na memoria
// TODO: tratar a seed nos nós
// TODO: Balancear a carga de trabalho entre os nós
// TODO: Usar primitivas coletivas do MPI
// TODO: Usar parallel/for/tasks/SIMD do OpenMP

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <N> <SEED> <T>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int N, SEED, T, *x, *y, *z;

    N = atoi(argv[1]); // Tamanho da matriz
    SEED = atoi(argv[2]); // Semente para gerar os valores aleatórios
    T = atoi(argv[3]); // Número de threads

    if (N <= 0 || SEED <= 0 || T <= 0) {
        printf("Invalid arguments, N and SEED must be positive numbers\n");
        return EXIT_FAILURE;
    }

    srand(SEED);

    omp_set_num_threads(T);

    int process_number, myrank;
    int src = 0, dest, message_tag;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &process_number);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Executa o processo 0
    if (myrank == 0) {
        // TODO: tratar o caso de N n caber na memória de apenas uma máquina
        x = createMatrix(N, N, SEED);
        y = createMatrix(N, N, SEED);
        z = createMatrix(N, N, SEED);

        // Envia as matrizes x, y, z para os demais processos
        MPI_Bcast(x, N*N, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(y, N*N, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(z, N*N, MPI_INT, src, MPI_COMM_WORLD);


        // printf("Distância de Manhattan mínima: %d (soma min: %d) e máxima: %d (soma max: %d).\n", min_manhattan, sum_min_manhattan, max_manhattan, sum_max_manhattan);
        // printf("Distância Euclidiana mínima: %.2lf (soma min: %.2lf) e máxima: %.2lf (soma max: %.2lf).\n", min_euclidean, sum_min_euclidean, max_euclidean, sum_max_euclidean);
  
    } else { // Executa os demais processos
    
        // Recebe as matrizes x, y, z do processo 0
        MPI_Bcast(x, N*N, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(y, N*N, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(z, N*N, MPI_INT, src, MPI_COMM_WORLD);

    }

    free(x);
    free(y);
    free(z);

    if (MPI_Finalize() != MPI_SUCCESS) {
        printf("Error on MPI_Finalize()\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// Cria uma matriz de tamanho rowXcolumn com valores aleatórios
int *createMatrix(int row, int column, int SEED) {
    int *matrix = (int *) malloc(row * column * sizeof(int*));

    for (int i = 0; i < row; i++)
        for (int j = 0; j < column; j++)
            matrix[i*column + j] = rand() % MAX_VALUE;

    return matrix;
}

void printMatrix(int *matrix, int row, int column) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++)
            printf("%d ", matrix[i*column + j]);

        printf("\n");
    }

    printf("\n");
}

// Libera a memória alocada para a matriz
void freeMatrix(int *matrix, int row) {
    free(matrix);
}