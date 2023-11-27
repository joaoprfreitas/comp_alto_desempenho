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
#include <string.h>

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
// TODO: Usar tasks
// TODO: tratar overflow nas somas
// TODO: processos e threads estão com cargas balanceadas?
// TODO: arrumar comandos de compilacao e add hosts

// Trata o caso de N não ser divisível inteiramente por n_processes
// Os processos com rank < NODE_MOD recebem um bloco de tamanho NODE_SIZE + 1
// 0 1 2 3
// [x, z, k]
// [p, z, _, _]
int getBlockSize(int myrank, int NODE_SIZE, int NODE_MOD) {
    return (myrank < NODE_MOD) ? NODE_SIZE + 1 : NODE_SIZE;
}

// 7 6 6 6

// myrank 0, n_processes 4, blockSize 7
// div = 7/4 = 1
// mod = 7%4 = 3
// 0 < 3 -> V -> div + 1 = 2
// 1 < 3 -> V -> div + 1 = 2
// 2 < 3 -> V -> div + 1 = 2
// 3 < 3 -> F -> div = 1

// myrank 1, n_processes 4, blockSize 6
// div = 6/4 = 1
// mod = 6%4 = 2
// 0 < 2 -> V -> div + 1 = 2
// 1 < 2 -> V -> div + 1 = 2
// 2 < 2 -> F -> div = 1
// 3 < 2 -> F -> div = 1


// myrank, n_processes, REAL_NODE_SIZE
int getTamanhoBloquinho(int myrank, int n_processes, int blockSize) {
    int div = blockSize / n_processes; 
    int mod = blockSize % n_processes; 

    if (div == 0) return 1; // Bloco mínimo é 1

    if (myrank < mod) { 
        return div + 1;
    } else {
        return div;
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
    // int src = 0;

    // Inicializa o MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int hostLenght;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(hostname, &hostLenght);

    // printf("Process %d of %d running on %s\n", myrank, n_processes, hostname);

    const int MATRIX_LENGTH = N * N;
    const int NODE_SIZE = (MATRIX_LENGTH / n_processes);
    const int NODE_MOD = (MATRIX_LENGTH % n_processes);

    const int REAL_NODE_SIZE = getBlockSize(myrank, NODE_SIZE, NODE_MOD);
    
    int *x = (int *) malloc(REAL_NODE_SIZE * sizeof(int));
    int *y = (int *) malloc(REAL_NODE_SIZE * sizeof(int));
    int *z = (int *) malloc(REAL_NODE_SIZE * sizeof(int));


    // * DAQUI PRA CIMA TA CERTO

    // TODO: tem q tratar dnv o +1 nesse divisao?
    int *aux = (int*) malloc((NODE_SIZE/n_processes + 1) * sizeof(int)); 

    int *x_recebido = (int *) malloc(NODE_SIZE * sizeof(int));
    int *y_recebido = (int *) malloc(NODE_SIZE * sizeof(int));
    int *z_recebido = (int *) malloc(NODE_SIZE * sizeof(int));

    int *dist_min_parcial_manhattan = (int *) malloc(NODE_SIZE * sizeof(int));
    int *dist_max_parcial_manhattan = (int *) malloc(NODE_SIZE * sizeof(int));
    double *dist_min_parcial_euclidiana = (double *) malloc(NODE_SIZE * sizeof(double));
    double *dist_max_parcial_euclidiana = (double *) malloc(NODE_SIZE * sizeof(double));

    int minManhattanGlobal = INT_MAX, maxManhattanGlobal = INT_MIN;
    double minEuclideanGlobal = INT_MAX, maxEuclideanGlobal = INT_MIN;
    double somaMinEuclidiana = 0;
    double somaMaxEuclidiana = 0;
    int somaMinManhattan = 0;
    int somaMaxManhattan = 0;

    for (int i = 0; i < NODE_SIZE; i++) {
        dist_min_parcial_manhattan[i] = INT_MAX;
        dist_max_parcial_manhattan[i] = INT_MIN;
        dist_min_parcial_euclidiana[i] = DBL_MAX;
        dist_max_parcial_euclidiana[i] = DBL_MIN;
    }

    if (myrank == MASTER_PROCESS) {
        int index_x0 = 0;
        int index_y0 = 0;
        int index_z0 = 0;

        int *x0 = (int*) malloc((REAL_NODE_SIZE) * sizeof(int));
        int *y0 = (int*) malloc((REAL_NODE_SIZE) * sizeof(int));
        int *z0 = (int*) malloc((REAL_NODE_SIZE) * sizeof(int));

        // printf("ENVIADO\n");

        int processo_atual = 0;
        int comecodoBloco = 0;

        for(int k = 0; k < n_processes; k++) {
            int tamBlocoEnviado = getBlockSize(k, NODE_SIZE, NODE_MOD);
            // printf("Bloco k = %d, tamBlocoEnviado = %d\n", k, tamBlocoEnviado);
            // printf("Valores enviados: ");

            // printf("Gerou %d valores para o bloco %d\n", tamBlocoEnviado, k);
            for (int i = 0; i < tamBlocoEnviado; i++) {
                x[i] = rand() % MAX_VALUE;
                // printf("x[%d] = %d ", i, x[i]);
            }
            // printf("\n");
            
            int comecoDoEnvio = comecodoBloco;
            int final = comecodoBloco + tamBlocoEnviado;
            int contagemDeProcessos = 0;
            
            while (comecoDoEnvio < final && contagemDeProcessos < n_processes) {
                int tamanhoBloquinho = 0;
                
                for (int j = comecoDoEnvio; j < final; j += n_processes) {
                    aux[(j - comecodoBloco) / n_processes] = x[j - comecodoBloco];
                    tamanhoBloquinho++;
                    // for (int z = 0; z < REAL_NODE_SIZE; z++) {
                    //     printf("%d ", aux[z]);
                    // }
                    // printf("\n");
                    // printf("Bloco k = %d, envia para processo %d, tamBloquinho = %d, j = %d, valor = %d\n", k, processo_atual, tamanhoBloquinho, j, aux[(j - comecodoBloco) / n_processes]);
                }

                comecoDoEnvio++;
                contagemDeProcessos++;
                
                if (processo_atual != 0) {
                    // printf("Envia para Processo %d - Tam Bloquinho %d - Aux: ", tamanhoBloquinho, processo_atual);
                    // for (int i = 0; i < REAL_NODE_SIZE; i++) {
                    //     printf("%d ", aux[i]);
                    // }
                    // printf("\n");
                    
                    printf("Enviando para processo %d: ", processo_atual);
                    for (int z = 0; z < tamanhoBloquinho; z++) {
                        printf("%d ", aux[z]);
                    }
                    printf("\n");

                    // printf("tam bloq enviado %d para o processo %d\n", tamanhoBloquinho, processo_atual);

                    MPI_Send(&tamanhoBloquinho, 1, MPI_INT, processo_atual, messageTag, MPI_COMM_WORLD);
                    messageTag++;
                    MPI_Send(aux, tamanhoBloquinho, MPI_INT, processo_atual, messageTag, MPI_COMM_WORLD);
                    messageTag--;
                    if (processo_atual == n_processes-1) {
                        messageTag+=2;
                    }        
                } else {
                    // TODO
                }
                
                processo_atual++;
                processo_atual %= n_processes;
            }

            comecodoBloco = final;
        }
        
        
        // messageTag = 0;
        for(int k = 0; k < n_processes; k++){
            for (int i = 0; i < NODE_SIZE; i++) {
                y[i] = rand() % MAX_VALUE;
            }
            for (int numero_processo = 1; numero_processo < n_processes; numero_processo++) {
                    for (int j = 0; j < NODE_SIZE/n_processes; j++) {
                        aux[j] = y[j * n_processes + numero_processo];
                    }
            
                    MPI_Send(aux, NODE_SIZE/n_processes, MPI_INT, numero_processo, messageTag, MPI_COMM_WORLD);   
            }
            messageTag++;
            for (int j = 0; j < NODE_SIZE/n_processes; j++) {
                y0[index_y0++] = y[j * n_processes];
            }
        }

        // messageTag = 0;
        
        for(int k = 0; k < n_processes; k++){
            for (int i = 0; i < NODE_SIZE; i++) {
                z[i] = rand() % MAX_VALUE;
            }
            for (int numero_processo = 1; numero_processo < n_processes; numero_processo++) {
                for (int j = 0; j < NODE_SIZE/n_processes; j++) {
                    aux[j] = z[j * n_processes + numero_processo];
                }
                MPI_Send(aux, NODE_SIZE/n_processes, MPI_INT, numero_processo, messageTag, MPI_COMM_WORLD);
            }
            messageTag++;
            for (int j = 0; j < NODE_SIZE/n_processes; j++) {
                z0[index_z0++] = z[j * n_processes];
            }
        }
  
        memcpy(x, x0, NODE_SIZE * sizeof *x);
        memcpy(y, y0, NODE_SIZE * sizeof *y);
        memcpy(z, z0, NODE_SIZE * sizeof *z);

        free(x0);
        free(y0);
        free(z0);
    } else {
        // TODO: o tamanho nodesize no receive ta certo?
        int mensageTag = 0;
        // int tamBlocoRecebido = getBlockSize(myrank, REAL_NODE_SIZE/n_processes, REAL_NODE_SIZE % n_processes);
        // int tamBlocoRecebido = getTamanhoBloquinho(myrank, n_processes, REAL_NODE_SIZE);
       

        // printf("Myrank: %d Tamanho bloquinho: %d, RealNodeSize = %d\n", myrank, tamBlocoRecebido, REAL_NODE_SIZE);


        int pos = 0;
        for (int i = 0; i < n_processes; i++) {
            int tamBloquinho = 0;
            MPI_Recv(&tamBloquinho, 1, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
            

            printf("Tam bloquinho %d recebido no rank %d\n", tamBloquinho, myrank);
            MPI_Recv(aux + pos, tamBloquinho, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
            
            pos += tamBloquinho;
        }

        for (int j = 0; j < REAL_NODE_SIZE; j++) {
            printf("%d ", aux[j]);
        }
        printf("\n");

        // int pos = 0;
        // for (int i = 0; i < n_processes; i++) {
        //     int tamBloquinho = 0;
        //     MPI_Recv(&tamBloquinho, 1, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
        //     printf("Tam bloquinho %d recebido no rank %d\n", tamBloquinho, myrank);
        //     int *teste = (int*) malloc(tamBloquinho * sizeof(int));
            
        //     for (int z = 0; z < tamBloquinho; z++) {
        //         teste[z] = -1;
        //     }

        //     // printf("Tamanho bloquinho: %d\n", tamBloquinho);
        //     // MPI_Recv(aux + pos, tamBloquinho, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
        //     MPI_Recv(teste, tamBloquinho, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
            
        //     // printf("Processo %d recebeu: ", myrank);
        //     for (int j = 0; j < tamBloquinho; j++) {
        //         printf("%d ", teste[j]);
        //     }
        //     printf("\n");
        //     // pos += tamBloquinho;
        //     free(teste);
        // }

        // int posicao = 0;
        // for (int i = 0; i < n_processes; i++) {
        //     if (tamBlocoRecebido == 1) {
        //         int valor = -1;
        //         MPI_Recv(&valor, 1, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);

        //         if (valor != -1) {
        //             x[posicao++] = valor;
        //             // printf("Valor recebido: %d\n", valor);
        //         }
        //     } else {
        //         MPI_Recv(x + posicao, tamBlocoRecebido, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
        //         posicao += tamBlocoRecebido;
        //     }
        // }

        // for (int i = 0; i < n_processes; i++) {
        //     if (myrank == 1) {
        //         // printf("Tag msg recebida: %d\n", mensageTag);
        //     }
        //     MPI_Recv(x + i * n_processes, NODE_SIZE, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
        // }

        // printf("Processo %d -- ", myrank);
        // for (int i = 0; i < REAL_NODE_SIZE; i++) {
        //     printf("X[%d] = %d ", i, x[i]);
        // }
        // printf("\n");
        
        for (int i = 0; i < n_processes; i++) {
            MPI_Recv(y+ i * n_processes, NODE_SIZE, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
        }

        for (int i = 0; i < n_processes; i++) {
            MPI_Recv(z + i * n_processes, NODE_SIZE, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
        }

    }
    
    //TODO: n deveria ser pares?
    MPI_Bcast(&x, NODE_SIZE, MPI_INT, myrank, MPI_COMM_WORLD);
    MPI_Bcast(&y, NODE_SIZE, MPI_INT, myrank, MPI_COMM_WORLD);
    MPI_Bcast(&z, NODE_SIZE, MPI_INT, myrank, MPI_COMM_WORLD);

    for (int i = 0; i < n_processes; i++) {
        // TODO: reduzir uma comunicacao nao recebendo de si mesmo
        // if (myrank != i) {
            MPI_Bcast(&x_recebido, NODE_SIZE, MPI_INT, i, MPI_COMM_WORLD);
            MPI_Bcast(&y_recebido, NODE_SIZE, MPI_INT, i, MPI_COMM_WORLD);
            MPI_Bcast(&z_recebido, NODE_SIZE, MPI_INT, i, MPI_COMM_WORLD);
            
            for (int j = 0; j < NODE_SIZE; j++) { // posição j -> menor distancia em relação ao vetor que tá no processo i
                // int minManhattanLocal_parcial = INT_MAX, maxManhattanLocal_parcial = INT_MIN;
                // double minEuclideanLocal_parcial = INT_MAX, maxEuclideanLocal_parcial = INT_MIN;

                // #pragma
                for (int k = 0; k < NODE_SIZE; k++) {
                    if ((j * n_processes + myrank) < (k * n_processes + i)) {
                        // K: dados do processo que ta chegando por broadcast
                        // J: dados do processo atual
                        int manhattanDistance = manhattan(x[j], y[j], z[j], x_recebido[k], y_recebido[k], z_recebido[k]);
                        double euclideanDistance = euclidean(x[j], y[j], z[j], x_recebido[k], y_recebido[k], z_recebido[k]);

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
                } // end k
            } // end j
        // }
    } // end i

    // fazer reduce entre todos os processos, pois temos 4 vetores de distancias
    // 

    if (myrank == 0) {
        MPI_Reduce(&dist_min_parcial_manhattan, &minManhattanGlobal, NODE_SIZE, MPI_INT, MPI_MIN, MASTER_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(&dist_max_parcial_manhattan, &maxManhattanGlobal, NODE_SIZE, MPI_INT, MPI_MAX, MASTER_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(&dist_min_parcial_euclidiana, &minEuclideanGlobal, NODE_SIZE, MPI_DOUBLE, MPI_MIN, MASTER_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(&dist_max_parcial_euclidiana, &maxEuclideanGlobal, NODE_SIZE, MPI_DOUBLE, MPI_MAX, MASTER_PROCESS, MPI_COMM_WORLD);

        MPI_Reduce(&dist_min_parcial_manhattan, &somaMinManhattan, NODE_SIZE, MPI_INT, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(&dist_max_parcial_manhattan, &somaMaxManhattan, NODE_SIZE, MPI_INT, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(&dist_min_parcial_euclidiana, &somaMinEuclidiana, NODE_SIZE, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(&dist_max_parcial_euclidiana, &somaMaxEuclidiana, NODE_SIZE, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);


        printf("Distância de Manhattan mínima: %d (soma min: %d) e máxima: %d (soma max: %d).\n", minManhattanGlobal, somaMinManhattan, maxManhattanGlobal, somaMaxManhattan);
        printf("Distância Euclidiana mínima: %.2lf (soma min: %.2lf) e máxima: %.2lf (soma max: %.2lf).\n", minEuclideanGlobal, somaMinEuclidiana, maxEuclideanGlobal, somaMaxEuclidiana);
    }

    free(x);
    free(y);
    free(z);
    free(x_recebido);
    free(y_recebido);
    free(z_recebido);
    free(aux);
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

void printMatrix(int *matrix, int row, int column) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++)
            printf("%d ", matrix[i*column + j]);

        printf("\n");
    }

    printf("\n");
}