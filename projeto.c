/*
 * Integrantes:
 * João Pedro Rodrigues Freitas - 11316552
 * Gabriel Akio Urakawa - 11795912
 * Renato Tadeu Theodoro Junior - 11796750
 * Samuel Victorio Bernacci - 12703455
 * 
 * Comandos para compilar:
 * mpicc main.c -o main -fopenmp -lm -Wall
 * mpirun -np 32 --hostfile hosts.txt ./main N SEED T
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#define MAX_VALUE 100
#define MASTER_PROCESS 0

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

// Envia a matriz de dados geradas pela seed para o processo correto
void generateMatrix(int MATRIX_LENGTH, int n_processes, int *matrix, int messageTag){
    int valorGerado;
    for (int i = 0; i < MATRIX_LENGTH; i++) {
        valorGerado = rand() % MAX_VALUE;
        if (i % n_processes == 0) {
            matrix[i / n_processes] = valorGerado;
        } else {
            //O envio é feito ponto a ponto, ou seja geramos um valor e enviamos para um processo.
            // É feito um envio circular para balancear os dados.
            MPI_Send(&valorGerado, 1, MPI_INT, i % n_processes, messageTag, MPI_COMM_WORLD);
        }
    }
}

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

    srand(SEED);

    omp_set_num_threads(T); // Define o número de threads
    omp_set_nested(1); // Permite o uso de threads dentro de threads

    int n_processes, myrank, messageTag = 0;

    // Inicializa o MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    const int MATRIX_LENGTH = N * N;
    const int NODE_SIZE = MATRIX_LENGTH / n_processes;
    const int NODE_MOD = MATRIX_LENGTH % n_processes;
    
    // Numero de elementos que cada processo vai receber
    // Se o numero de elementos nao for divisivel pelo numero de processos, os primeiros processos
    // recebem um elemento a mais cada um até esgotar o NODE_MOD
    const int REAL_NODE_SIZE = NODE_SIZE + (myrank < NODE_MOD ? 1 : 0);
    
    int *x = (int *) malloc(REAL_NODE_SIZE * sizeof(int));
    int *y = (int *) malloc(REAL_NODE_SIZE * sizeof(int));
    int *z = (int *) malloc(REAL_NODE_SIZE * sizeof(int));

    // Vetores que armazenam as distâncias mínimas e máximas parciais de cada processo
    int *dist_min_parcial_manhattan = (int *) malloc((NODE_SIZE + 1) * sizeof(int));
    int *dist_max_parcial_manhattan = (int *) malloc((NODE_SIZE + 1) * sizeof(int));
    double *dist_min_parcial_euclidiana = (double *) malloc((NODE_SIZE + 1) * sizeof(double));
    double *dist_max_parcial_euclidiana = (double *) malloc((NODE_SIZE + 1) * sizeof(double));

    // Inicializa os vetores parciais com valores extremos
    for (int i = 0; i < NODE_SIZE + 1; i++) {
        dist_min_parcial_manhattan[i] = INT_MAX;
        dist_max_parcial_manhattan[i] = INT_MIN;
        dist_min_parcial_euclidiana[i] = DBL_MAX;
        dist_max_parcial_euclidiana[i] = DBL_MIN;
    }
    
    // Inicializa os valores globais para cada processo
    int minManhattanGlobal = INT_MAX, maxManhattanGlobal = INT_MIN;
    double minEuclideanGlobal = INT_MAX, maxEuclideanGlobal = INT_MIN;
    double somaMinEuclidiana = 0;
    double somaMaxEuclidiana = 0;
    int somaMinManhattan = 0;
    int somaMaxManhattan = 0;

    if (myrank == MASTER_PROCESS) {
        // O processo mestre gera os valores e os distribui entre todos os processos
        // Enquanto mantém uma parcela dos valores para ser processada pelo processo 0 também
        generateMatrix(MATRIX_LENGTH, n_processes, x, messageTag);
        generateMatrix(MATRIX_LENGTH, n_processes, y, messageTag);
        generateMatrix(MATRIX_LENGTH, n_processes, z, messageTag);
    } else {
        int x_recebido, y_recebido, z_recebido;
        // Recebimento e armazenamento dos dados nos diferentes nós
        for (int i = 0; i < REAL_NODE_SIZE; i++) {
            MPI_Recv(&x_recebido, 1, MPI_INT, MASTER_PROCESS, messageTag, MPI_COMM_WORLD, &status);
            x[i] = x_recebido;
        }
    
        for (int i = 0; i < REAL_NODE_SIZE; i++) {
            MPI_Recv(&y_recebido, 1, MPI_INT, MASTER_PROCESS, messageTag, MPI_COMM_WORLD, &status);
            y[i] = y_recebido;
        }
        
        for (int i = 0; i < REAL_NODE_SIZE; i++) {
            MPI_Recv(&z_recebido, 1, MPI_INT, MASTER_PROCESS, messageTag, MPI_COMM_WORLD, &status);
            z[i] = z_recebido;
        }
    }

    // Cálculo das distâncias e compartilhamento de dados 
    for (int i = 0; i < n_processes; i++) {
        int NODE_SIZE_RECEIVED = NODE_SIZE + (i < NODE_MOD ? 1 : 0);

        // armazenamento temporário
        int *x_recebido = (int *) malloc(NODE_SIZE_RECEIVED * sizeof(int));
        int *y_recebido = (int *) malloc(NODE_SIZE_RECEIVED * sizeof(int));
        int *z_recebido = (int *) malloc(NODE_SIZE_RECEIVED * sizeof(int));

        // Envio por meio de broadcasts dos dados de um nó para todos os outros para cálculo das distâncias
        if (myrank == i) {
            MPI_Bcast(x, REAL_NODE_SIZE, MPI_INT, myrank, MPI_COMM_WORLD);
            MPI_Bcast(y, REAL_NODE_SIZE, MPI_INT, myrank, MPI_COMM_WORLD);
            MPI_Bcast(z, REAL_NODE_SIZE, MPI_INT, myrank, MPI_COMM_WORLD);
            // Tratamento do caso onde o brodcast é "recebido" pelo nó que o enviou
            memcpy(x_recebido, x, REAL_NODE_SIZE * sizeof(int));
            memcpy(y_recebido, y, REAL_NODE_SIZE * sizeof(int));
            memcpy(z_recebido, z, REAL_NODE_SIZE * sizeof(int));
        } 
        // Recebimento dos broadcasts
        else {
            MPI_Bcast(x_recebido, NODE_SIZE_RECEIVED, MPI_INT, i, MPI_COMM_WORLD);
            MPI_Bcast(y_recebido, NODE_SIZE_RECEIVED, MPI_INT, i, MPI_COMM_WORLD);
            MPI_Bcast(z_recebido, NODE_SIZE_RECEIVED, MPI_INT, i, MPI_COMM_WORLD);
        }
        
        // calculo das distâncias
        for (int j = 0; j < REAL_NODE_SIZE; j++) { 
            // Claúsula SIMD é utilizada no cálculo das distâncias, facilitando o processamento de extensivas operações iguais
            #pragma omp simd
            for (int k = 0; k < NODE_SIZE_RECEIVED; k++) {  
                if ((j * n_processes + myrank) < (k * n_processes + i)) {
                    int manhattanDistance = manhattan(x[j], y[j], z[j], x_recebido[k], y_recebido[k], z_recebido[k]);
                    double euclideanDistance = euclidean(x[j], y[j], z[j], x_recebido[k], y_recebido[k], z_recebido[k]);
                    
                    // Atualização dos mínimos e máximos parciais(referentes apenas aos dados do bloco que está no nó)
                    if (manhattanDistance < dist_min_parcial_manhattan[j]) {
                        dist_min_parcial_manhattan[j] = manhattanDistance;
                    }

                    if (manhattanDistance > dist_max_parcial_manhattan[j]) {
                        dist_max_parcial_manhattan[j] = manhattanDistance;
                    }

                    if (euclideanDistance < dist_min_parcial_euclidiana[j]) {
                        dist_min_parcial_euclidiana[j] = euclideanDistance;
                    }

                    if (euclideanDistance > dist_max_parcial_euclidiana[j]) {
                        dist_max_parcial_euclidiana[j] = euclideanDistance;
                    }
                }
            } 
        } 
        
        free(x_recebido);
        free(y_recebido);
        free(z_recebido);
    }

    int minLocalAbsolutoManhattan = INT_MAX, maxLocalAbsolutoManhattan = INT_MIN, sumMinLocalManhattan = 0, sumMaxLocalManhattan = 0; 
    double minLocalAbsolutoEuclidiana = INT_MAX, maxLocalAbsolutoEuclidiana = INT_MIN, sumMinLocalEuclidiana = 0, sumMaxLocalEuclidiana =0;
    // Reduction é utilizada para condensar os vários mínimos e máximos parciais para as suas versões absolutas, porém ainda locais
    #pragma omp parallel for \
        reduction(min: minLocalAbsolutoManhattan, minLocalAbsolutoEuclidiana) \
        reduction(max: maxLocalAbsolutoManhattan, maxLocalAbsolutoEuclidiana) \
        reduction(+: sumMinLocalManhattan, sumMaxLocalManhattan, sumMinLocalEuclidiana, sumMaxLocalEuclidiana) 
    for (int i = 0; i < REAL_NODE_SIZE; i++) {
        
        // Tasks são utilizadas para incrementar o paralelismo, permitindo a comparação simultânea de  diversas váriais diferentes,
        // Sem relação umas com as outras
        #pragma omp task shared(minLocalAbsolutoManhattan)
        {
            if (dist_min_parcial_manhattan[i] < minLocalAbsolutoManhattan) {
                minLocalAbsolutoManhattan = dist_min_parcial_manhattan[i];
            }
        }
        
        #pragma omp task shared(maxLocalAbsolutoManhattan)
        {
            if (dist_max_parcial_manhattan[i] > maxLocalAbsolutoManhattan) {
                maxLocalAbsolutoManhattan = dist_max_parcial_manhattan[i];
            }
        }

        #pragma omp task shared(minLocalAbsolutoEuclidiana)
        {
            if (dist_min_parcial_euclidiana[i] < minLocalAbsolutoEuclidiana) {
                minLocalAbsolutoEuclidiana = dist_min_parcial_euclidiana[i];
            }
        }

        #pragma omp task shared(maxLocalAbsolutoEuclidiana)
        {
            if (dist_max_parcial_euclidiana[i] > maxLocalAbsolutoEuclidiana) {
                maxLocalAbsolutoEuclidiana = dist_max_parcial_euclidiana[i] ;
            }
        }

        // Taskwait é utilizado para evitar race condition entre tasks de diferentes iterações
        #pragma omp taskwait
        {
            // atualiza a soma dos valores máximos e mínimos
            if (dist_min_parcial_manhattan[i] != INT_MAX) {
                sumMaxLocalManhattan += dist_max_parcial_manhattan[i];
                sumMinLocalEuclidiana += dist_min_parcial_euclidiana[i];
                sumMinLocalManhattan += dist_min_parcial_manhattan[i];
                sumMaxLocalEuclidiana += dist_max_parcial_euclidiana[i];
            }
        }

    }

    // fazer reduce entre todos os processos, pois temos 4 vetores de distancias 
    // Reductions utilizados para finalmente condensar os mínimos, máximos e somas locais absolutas em suas versões globais

    MPI_Reduce(&minLocalAbsolutoManhattan, &minManhattanGlobal, 1, MPI_INT, MPI_MIN, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&maxLocalAbsolutoManhattan, &maxManhattanGlobal, 1, MPI_INT, MPI_MAX, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&minLocalAbsolutoEuclidiana, &minEuclideanGlobal, 1, MPI_DOUBLE, MPI_MIN, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&maxLocalAbsolutoEuclidiana, &maxEuclideanGlobal, 1, MPI_DOUBLE, MPI_MAX, MASTER_PROCESS, MPI_COMM_WORLD);

    MPI_Reduce(&sumMinLocalManhattan, &somaMinManhattan, 1, MPI_INT, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&sumMaxLocalManhattan, &somaMaxManhattan, 1, MPI_INT, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&sumMinLocalEuclidiana, &somaMinEuclidiana, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&sumMaxLocalEuclidiana, &somaMaxEuclidiana, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);

    // Prints dos resultados finais
    if (myrank == MASTER_PROCESS) {
        printf("Distância de Manhattan mínima: %d (soma min: %d) e máxima: %d (soma max: %d).\n", minManhattanGlobal, somaMinManhattan, maxManhattanGlobal, somaMaxManhattan);
        printf("Distância Euclidiana mínima: %.2lf (soma min: %.2lf) e máxima: %.2lf (soma max: %.2lf).\n", minEuclideanGlobal, somaMinEuclidiana, maxEuclideanGlobal, somaMaxEuclidiana);
    }

    free(x);
    free(y);
    free(z);
    free(dist_min_parcial_manhattan);
    free(dist_max_parcial_manhattan);
    free(dist_min_parcial_euclidiana);
    free(dist_max_parcial_euclidiana);

    // Encerra o MPI
    if (MPI_Finalize() != MPI_SUCCESS) {
        printf("Error on MPI_Finalize(), process: %d\n", myrank);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

