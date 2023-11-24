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
#include <float.h>
#include <limits.h>

#define MAX_VALUE 100

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

    int *x = (int *) malloc(N * N * sizeof(int*));
    int *y = (int *) malloc(N * N * sizeof(int*));
    int *z = (int *) malloc(N * N * sizeof(int*));

    srand(SEED);

    omp_set_num_threads(T); // Define o número de threads
    omp_set_nested(1); // Permite o uso de threads dentro de threads

    int n_processes, myrank, src = 0, messageTag = 0;

    // Inicializa o MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int totalMinManhattan = 0, totalMaxManhattan = 0;
    int minGlobalManhattan = INT_MAX, maxGlobalManhattan = INT_MIN;
    double totalMinEuclidean = 0.0, totalMaxEuclidean = 0.0;
    double minGlobalEuclidean = DBL_MAX, maxGlobalEuclidean = DBL_MIN;

    int *manhattanInfos = (int *) malloc(4 * sizeof(int)); // min, max, sum_min, sum_max
    double *euclideanInfos = (double *) malloc(4 * sizeof(double)); // min, max, sum_min, sum_max

    int *receiveBufferManhattan = (int *) malloc((n_processes - 1) * sizeof(int)); // min, max, sum_min, sum_max
    double *receiveBufferEuclidean = (double *) malloc((n_processes - 1) * sizeof(double)); // min, max, sum_min, sum_max

    int receiveBufferSize = (n_processes - 1) * 4;

    int *startEnd = (int *) malloc(2 * sizeof(int)); // start, end

    // Executa o processo 0
    if (myrank == 0) { // TODO: da pra colocar o processo 0 para executar também
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
        MPI_Bcast(x, N*N, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(y, N*N, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(z, N*N, MPI_INT, src, MPI_COMM_WORLD);

        // TODO: o send deve ser bloqueante
        // TODO: lógica de balanceamento de carga

        int chunk = N * N / (n_processes - 1);
        int currentStart = 0, currentEnd = 0;
        for (int i = 1; i < n_processes; i++) {
            // send inicio e fim do for

            if (i == n_processes - 1) {
                currentEnd = N * N;
            } else {
                currentEnd += chunk;
            }

            startEnd[0] = currentStart;
            startEnd[1] = currentEnd;

            MPI_Send(startEnd, 2, MPI_INT, i, messageTag, MPI_COMM_WORLD);

            currentStart = startEnd[1];
        }

        // TODO: receber os dados dos outros processos
        MPI_Gather(manhattanInfos, 4, MPI_INT, receiveBufferManhattan, 4, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Gather(euclideanInfos, 4, MPI_DOUBLE, receiveBufferEuclidean, 4, MPI_DOUBLE, src, MPI_COMM_WORLD);

        // printf("Infos recebidas pelo processo 0\n");
        for (int i = 0; i < n_processes - 1; i++) {
            printf("RECEB Manhattan::: Min local: %d, max local: %d, soma min local: %d, soma max local: %d\n", receiveBufferManhattan[i * 4], receiveBufferManhattan[i * 4 + 1], receiveBufferManhattan[i * 4 + 2], receiveBufferManhattan[i * 4 + 3]);
            printf("RECEB Euclidean::: Min local: %.2lf, max local: %.2lf, soma min local: %.2lf, soma max local: %.2lf\n", receiveBufferEuclidean[i * 4], receiveBufferEuclidean[i * 4 + 1], receiveBufferEuclidean[i * 4 + 2], receiveBufferEuclidean[i * 4 + 3]);
        }

        // TODO: colocar no omp???
        for (int i = 0; i < n_processes - 1; i += 4) {
            if (minGlobalManhattan > manhattanInfos[i]) {
                minGlobalManhattan = manhattanInfos[i];
            }

            if (maxGlobalManhattan < manhattanInfos[i + 1]) {
                maxGlobalManhattan = manhattanInfos[i + 1];
            }

            if (minGlobalEuclidean > euclideanInfos[i]) {
                minGlobalEuclidean = euclideanInfos[i];
            }

            if (maxGlobalEuclidean < euclideanInfos[i + 1]) {
                maxGlobalEuclidean = euclideanInfos[i + 1];
            }

            totalMinManhattan += manhattanInfos[i + 2];
            totalMaxManhattan += manhattanInfos[i + 3];
            totalMinEuclidean += euclideanInfos[i + 2];
            totalMaxEuclidean += euclideanInfos[i + 3];
        }

        // printf("Distância de Manhattan mínima: %d (soma min: %d) e máxima: %d (soma max: %d).\n", minGlobalManhattan, totalMinManhattan, maxGlobalManhattan, totalMaxManhattan);
        // printf("Distância Euclidiana mínima: %.2lf (soma min: %.2lf) e máxima: %.2lf (soma max: %.2lf).\n", minGlobalEuclidean, totalMinEuclidean, maxGlobalEuclidean, totalMaxEuclidean);
  
    } else { // Executa os demais processos
    
        // Recebe as matrizes x, y, z do processo 0
        MPI_Bcast(x, N*N, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(y, N*N, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Bcast(z, N*N, MPI_INT, src, MPI_COMM_WORLD);

        // TODO: o receive deve ser bloqueante
        MPI_Recv(startEnd, 2, MPI_INT, src, messageTag, MPI_COMM_WORLD, &status);

        // receive inicio e fim do for
        // int start, end;

        // obtém os maximos e minimos do ponto (x[i], y[i], z[i])
        // TODO: terminar o for externo
        #pragma omp parallel for \
            reduction(min: minGlobalManhattan, minGlobalEuclidean) \
            reduction(max: maxGlobalManhattan, maxGlobalEuclidean) \
            reduction(+: totalMaxEuclidean, totalMaxManhattan, totalMinEuclidean, totalMinManhattan)
        for (int i = startEnd[0]; i < startEnd[1]; i++) {
            int minManhattanLocal = INT_MAX, maxManhattanLocal = INT_MIN;
            double minEuclideanLocal = INT_MAX, maxEuclideanLocal = INT_MIN;

            // obtém as distâncias de manhattan e euclidiana entre o ponto (x[i], y[i], z[i]) e todos os outros pontos
            #pragma omp parallel for reduction(min: minManhattanLocal, minEuclideanLocal) reduction(max: maxManhattanLocal, maxEuclideanLocal)
            for (int j = i + 1; j < N * N; j++) {
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
        
        manhattanInfos[0] = minGlobalManhattan;
        manhattanInfos[1] = maxGlobalManhattan;
        manhattanInfos[2] = totalMinManhattan;
        manhattanInfos[3] = totalMaxManhattan;

        euclideanInfos[0] = minGlobalEuclidean;
        euclideanInfos[1] = maxGlobalEuclidean;
        euclideanInfos[2] = totalMinEuclidean;
        euclideanInfos[3] = totalMaxEuclidean;

        printf("ENVIO Manhattan %d::: Min local: %d, max local: %d, soma min local: %d, soma max local: %d\n", myrank, minGlobalManhattan, maxGlobalManhattan, totalMinManhattan, totalMaxManhattan);
        printf("ENVIO Euclidean %d::: Min local: %.2lf, max local: %.2lf, soma min local: %.2lf, soma max local: %.2lf\n", myrank, minGlobalEuclidean, maxGlobalEuclidean, totalMinEuclidean, totalMaxEuclidean);
        
        // * ATÉ AQUI TA CERTO





        MPI_Gather(manhattanInfos, 4, MPI_INT, receiveBufferManhattan, 4, MPI_INT, src, MPI_COMM_WORLD);
        MPI_Gather(euclideanInfos, 4, MPI_DOUBLE, receiveBufferEuclidean, 4, MPI_DOUBLE, src, MPI_COMM_WORLD);
    } // end else

    free(x);
    free(y);
    free(z);
    free(manhattanInfos);
    free(euclideanInfos);
    free(startEnd);
    free(receiveBufferEuclidean);
    free(receiveBufferManhattan);

    // Encerra o MPI
    if (MPI_Finalize() != MPI_SUCCESS) {
        printf("Error on MPI_Finalize()\n");
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