/* Integrantes:
 * João Pedro Rodrigues Freitas - 11316552
 * Gabriel Akio Urakawa - 11795912
 * Renato Tadeu Theodoro Junior - 11796750
 * Samuel Victorio Bernacci - 12703455
 * 
 * Comandos para compilar:
 * mpicc main.c -o main -fopenmp -lm -Wall
 * mpirun --oversubscribe -np 4 ./main N SEED T
 * 
 * Caso prefira, deixamos um Makefile de teste
 * basta edita-lo com os valores desejados e executar:
 * make mpi
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#define MAX_VALUE 100
#define MASTER_PROCESS 0

// TODO: alterar os tipos para esses:
typedef unsigned long long int ulli;
typedef long double ld;

// Retorna a distância euclidiana entre dois pontos (x, y, z)
double euclidean(int x1, int y1, int z1, int x2, int y2, int z2) {
    int x = x1 - x2;
    int y = y1 - y2;
    int z = z1 - z2;
    return sqrt(x*x + y*y + z*z);
}

// Retorna a distância de manhattan entre dois pontos (x, y, z)
int manhattan(int x1, int y1, int z1, int x2, int y2, int z2) {
    return abs(x1 - x2) + abs(y1 - y2) + abs(z1 - z2);
}

void populateMatrix(int *matrix, int row, int column);
void printMatrix(int *matrix, int row, int column);

// TODO: tratar o caso de N n caber na memoria
// TODO: tratar a seed nos nós
// TODO: Balancear a carga de trabalho entre os nós
// TODO: Usar primitivas coletivas do MPI
// TODO: Usar parallel/for/tasks/SIMD do OpenMP
// TODO: tratar overflow nas somas

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <N> <SEED> <T>\n", argv[0]);
        return EXIT_FAILURE;
    }
    int N = atoi(argv[1]); // Tamanho da matriz
    int SEED = atoi(argv[2]); // Semente para gerar os valores aleatórios
    int T = atoi(argv[3]); // Número de threads

    if (N <= 0 || SEED <= 0 || T <= 0) {
        printf("Invalid arguments: N, SEED and T must be positive\n");
        return EXIT_FAILURE;
    }

    int *x = (int *) malloc(N * N * sizeof(int));
    int *y = (int *) malloc(N * N * sizeof(int));
    int *z = (int *) malloc(N * N * sizeof(int));

    const int MATRIX_SIZE = N * N;

    srand(SEED);

    omp_set_num_threads(T); // Define o número de threads
    omp_set_nested(1); // Permite o uso de threads dentro de threads

    int n_processes, myrank, src = 0, messageTag = 0;

    // Inicializa o MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int hostLenght;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(hostname, &hostLenght);

    printf("Process %d of %d running on %s\n", myrank, n_processes, hostname);

    int totalMinManhattan = 0, totalMaxManhattan = 0;
    int minGlobalManhattan = INT_MAX, maxGlobalManhattan = INT_MIN;
    double totalMinEuclidean = 0.0, totalMaxEuclidean = 0.0;
    double minGlobalEuclidean = DBL_MAX, maxGlobalEuclidean = DBL_MIN;

    // Buffer para armazenar os resultados parciais de cada processo
    int *manhattanInfos = (int *) malloc(4 * sizeof(int)); // min, max, sum_min, sum_max
    double *euclideanInfos = (double *) malloc(4 * sizeof(double)); // min, max, sum_min, sum_max

    // Tamanho do buffer para receber os resultados parciais dos processos
    int receiveBufferSize = (n_processes) * 4; // sao 4 valores por processo: min, max, sum_min, sum_max

    // Buffer para receber os resultados parciais dos processos
    int *receiveBufferManhattan = (int *) malloc((receiveBufferSize) * sizeof(int)); // min, max, sum_min, sum_max
    double *receiveBufferEuclidean = (double *) malloc((receiveBufferSize) * sizeof(double)); // min, max, sum_min, sum_max

    // Buffer para enviar a carga de trabalho para os processos [inicio, fim]
    int *startEnd = (int *) malloc(2 * sizeof(int));

    // Executa o processo 0
    if (myrank == MASTER_PROCESS) {
        populateMatrix(x, N, N);
        populateMatrix(y, N, N);
        populateMatrix(z, N, N);

        // TODO: tratar o caso de N n caber na memória de apenas uma máquina
        if (x == NULL || y == NULL || z == NULL) {
            // TODO
            printf("Error on malloc()\n");
            return EXIT_FAILURE;
        }

        // Envia as matrizes x, y, z para os demais processos
        MPI_Bcast(x, MATRIX_SIZE, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(y, MATRIX_SIZE, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(z, MATRIX_SIZE, MPI_INT, src, MPI_COMM_WORLD);

        // Balanceamento de carga
        // Divide a matriz em CHUNKS de tamanho igual e envia para os processos
        // para que cada um execute uma parte
        // O ultimo processo recebe CHUNK + resto da divisao
        int chunk = MATRIX_SIZE / (n_processes);
        int currentStart = chunk, currentEnd = chunk;
        for (int i = 1; i < n_processes; i++) {
            if (i == n_processes - 1) {
                currentEnd = MATRIX_SIZE;
            } else {
                currentEnd += chunk;
            }

            startEnd[0] = currentStart;
            startEnd[1] = currentEnd;
            
            MPI_Send(startEnd, 2, MPI_INT, i, messageTag, MPI_COMM_WORLD);

            currentStart = startEnd[1];
        }

        // O processo 0 também executa a parte dele, que vai de 0 a chunk
        startEnd[0] = 0;
        startEnd[1] = chunk;
    } else { // Demais processos
        // Recebe as matrizes x, y, z do processo 0
        MPI_Bcast(x, MATRIX_SIZE, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(y, MATRIX_SIZE, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(z, MATRIX_SIZE, MPI_INT, src, MPI_COMM_WORLD);

        // Recebe a sua carga de trabalho
        MPI_Recv(startEnd, 2, MPI_INT, src, messageTag, MPI_COMM_WORLD, &status);
    }

    // Todos os processos realizam operacoes em sua carga de trabalho

    // obtém os maximos e minimos do ponto (x[i], y[i], z[i])
    // TODO: terminar o for externo
    // TODO: colocar um schedule?
    #pragma omp parallel for \
        reduction(min: minGlobalManhattan, minGlobalEuclidean) \
        reduction(max: maxGlobalManhattan, maxGlobalEuclidean) \
        reduction(+: totalMaxEuclidean, totalMaxManhattan, totalMinEuclidean, totalMinManhattan)
    for (int i = startEnd[0]; i < startEnd[1]; i++) {
        int minManhattanLocal = INT_MAX, maxManhattanLocal = INT_MIN;
        double minEuclideanLocal = INT_MAX, maxEuclideanLocal = INT_MIN;

        // obtém as distâncias de manhattan e euclidiana entre o ponto (x[i], y[i], z[i]) e todos os outros pontos
        #pragma omp parallel for reduction(min: minManhattanLocal, minEuclideanLocal) reduction(max: maxManhattanLocal, maxEuclideanLocal)
        for (int j = i + 1; j < MATRIX_SIZE; j++) {
            // manhattan_distance e euclidean_distance são variáveis locais de cada thread
            int manhattanDistance = manhattan(x[i], y[i], z[i], x[j], y[j], z[j]);
            double euclideanDistance = euclidean(x[i], y[i], z[i], x[j], y[j], z[j]);

            // atualiza os valores máximos e mínimos
            if (manhattanDistance < minManhattanLocal) {
                minManhattanLocal = manhattanDistance;
            }

            if (manhattanDistance > maxManhattanLocal) {
                maxManhattanLocal = manhattanDistance;
            }

            if (euclideanDistance < minEuclideanLocal) {
                minEuclideanLocal = euclideanDistance;
            }

            if (euclideanDistance > maxEuclideanLocal) {
                maxEuclideanLocal = euclideanDistance;
            }
        } // end for j

        // atualiza os valores máximos e mínimos globais
        if (minManhattanLocal < minGlobalManhattan) {
            minGlobalManhattan = minManhattanLocal;
        }

        if (maxManhattanLocal > maxGlobalManhattan) {
            maxGlobalManhattan = maxManhattanLocal;
        }

        if (minEuclideanLocal < minGlobalEuclidean) {
            minGlobalEuclidean = minEuclideanLocal;
        }

        if (maxEuclideanLocal > maxGlobalEuclidean) {
            maxGlobalEuclidean = maxEuclideanLocal;
        }

        // atualiza a soma dos valores máximos e mínimos
        if (minManhattanLocal != INT_MAX && minEuclideanLocal != DBL_MAX) {
            totalMaxManhattan += maxManhattanLocal;
            totalMaxEuclidean += maxEuclideanLocal;
            totalMinManhattan += minManhattanLocal;
            totalMinEuclidean += minEuclideanLocal;
        }
    } // end for i

    if (myrank == MASTER_PROCESS) {
        // Armazena no buffer os resultados do processo processo mestre, que também executou a parte dele
        receiveBufferManhattan[0] = minGlobalManhattan;
        receiveBufferManhattan[1] = maxGlobalManhattan;
        receiveBufferManhattan[2] = totalMinManhattan;
        receiveBufferManhattan[3] = totalMaxManhattan;

        receiveBufferEuclidean[0] = minGlobalEuclidean;
        receiveBufferEuclidean[1] = maxGlobalEuclidean;
        receiveBufferEuclidean[2] = totalMinEuclidean;
        receiveBufferEuclidean[3] = totalMaxEuclidean;

        // Recebe os dados parciais dos outros processos
        MPI_Gather(manhattanInfos, 4, MPI_INT, receiveBufferManhattan, 4, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Gather(euclideanInfos, 4, MPI_DOUBLE, receiveBufferEuclidean, 4, MPI_DOUBLE, src, MPI_COMM_WORLD);

        // TODO: colocar no omp???
        // Calcula os valores globais e soma os valores parciais
        for (int i = 4; i < receiveBufferSize; i += 4) {
            if (minGlobalManhattan > receiveBufferManhattan[i]) {
                minGlobalManhattan = receiveBufferManhattan[i];
            }

            if (maxGlobalManhattan < receiveBufferManhattan[i + 1]) {
                maxGlobalManhattan = receiveBufferManhattan[i + 1];
            }

            if (minGlobalEuclidean > receiveBufferEuclidean[i]) {
                minGlobalEuclidean = receiveBufferEuclidean[i];
            }

            if (maxGlobalEuclidean < receiveBufferEuclidean[i + 1]) {
                maxGlobalEuclidean = receiveBufferEuclidean[i + 1];
            }

            totalMinManhattan += receiveBufferManhattan[i + 2];
            totalMaxManhattan += receiveBufferManhattan[i + 3];
            totalMinEuclidean += receiveBufferEuclidean[i + 2];
            totalMaxEuclidean += receiveBufferEuclidean[i + 3];
        }

        printf("Distância de Manhattan mínima: %d (soma min: %d) e máxima: %d (soma max: %d).\n", minGlobalManhattan, totalMinManhattan, maxGlobalManhattan, totalMaxManhattan);
        printf("Distância Euclidiana mínima: %.2lf (soma min: %.2lf) e máxima: %.2lf (soma max: %.2lf).\n", minGlobalEuclidean, totalMinEuclidean, maxGlobalEuclidean, totalMaxEuclidean);
    } else {
        // Os demais processos enviam seus resultados parciais para o processo mestre

        // Armazena no buffer os resultados do processo
        manhattanInfos[0] = minGlobalManhattan;
        manhattanInfos[1] = maxGlobalManhattan;
        manhattanInfos[2] = totalMinManhattan;
        manhattanInfos[3] = totalMaxManhattan;

        euclideanInfos[0] = minGlobalEuclidean;
        euclideanInfos[1] = maxGlobalEuclidean;
        euclideanInfos[2] = totalMinEuclidean;
        euclideanInfos[3] = totalMaxEuclidean;

        // Envia os dados parciais para o processo mestre
        MPI_Gather(manhattanInfos, 4, MPI_INT, receiveBufferManhattan, 4, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Gather(euclideanInfos, 4, MPI_DOUBLE, receiveBufferEuclidean, 4, MPI_DOUBLE, src, MPI_COMM_WORLD);
    }

    free(x);
    free(y);
    free(z);
    free(startEnd);
    free(manhattanInfos);
    free(euclideanInfos);
    free(receiveBufferEuclidean);
    free(receiveBufferManhattan);
    fflush(0);

    // Encerra o MPI
    if (MPI_Finalize() != MPI_SUCCESS) {
        printf("Error on MPI_Finalize(), process: %d\n", myrank);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// Cria uma matriz de tamanho rowXcolumn com valores aleatórios
void populateMatrix(int *matrix, int row, int column) {
    for (int i = 0; i < row; i++)
        for (int j = 0; j < column; j++)
            matrix[i*column + j] = rand() % MAX_VALUE;
}

void printMatrix(int *matrix, int row, int column) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++)
            printf("%d ", matrix[i*column + j]);

        printf("\n");
    }

    printf("\n");
}