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

    // if(N % n_processes != 0){
    //     printf("Não funciona :(\n");
    //     if (MPI_Finalize() != MPI_SUCCESS) {
    //         printf("Error on MPI_Finalize(), process: %d\n", myrank);
    //         return EXIT_FAILURE;
    //     }
    // }


    int hostLenght;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(hostname, &hostLenght);

    // printf("Process %d of %d running on %s\n", myrank, n_processes, hostname);

    const int MATRIX_LENGTH = N * N;

    // adicionamos o +1 para contornar o arredondamento
    int NODE_SIZE = (MATRIX_LENGTH / n_processes);
    int NODE_SIZE_0 = NODE_SIZE + MATRIX_LENGTH % n_processes;
    printf("Node Size %d e Node_0 %d\n", NODE_SIZE, NODE_SIZE_0);
    int *x = (int *) malloc(NODE_SIZE_0 * sizeof(int));
    int *y = (int *) malloc(NODE_SIZE_0 * sizeof(int));
    int *z = (int *) malloc(NODE_SIZE_0 * sizeof(int));
    int *aux = (int*) malloc((NODE_SIZE_0/n_processes) * sizeof(int));

    
    int *dist_min_parcial_manhattan = (int *) malloc(NODE_SIZE_0 * sizeof(int));
    int *dist_max_parcial_manhattan = (int *) malloc(NODE_SIZE_0 * sizeof(int));
    double *dist_min_parcial_euclidiana = (double *) malloc(NODE_SIZE_0 * sizeof(double));
    double *dist_max_parcial_euclidiana = (double *) malloc(NODE_SIZE_0 * sizeof(double));

    int minManhattanGlobal = INT_MAX, maxManhattanGlobal = INT_MIN;
    double minEuclideanGlobal = INT_MAX, maxEuclideanGlobal = INT_MIN;
    double somaMinEuclidiana = 0;
    double somaMaxEuclidiana = 0;
    int somaMinManhattan = 0;
    int somaMaxManhattan = 0;

    for (int i = 0; i < NODE_SIZE_0; i++) {
        dist_min_parcial_manhattan[i] = INT_MAX;
        dist_max_parcial_manhattan[i] = INT_MIN;
        dist_min_parcial_euclidiana[i] = INT_MAX;
        dist_max_parcial_euclidiana[i] = INT_MIN;
    }

    if (myrank == MASTER_PROCESS) {
        int index_x0 = 0;
        int index_y0 = 0;
        int index_z0 = 0;
        int *x0 = (int*) malloc((NODE_SIZE_0) * sizeof(int));
        int *y0 = (int*) malloc((NODE_SIZE_0) * sizeof(int));
        int *z0 = (int*) malloc((NODE_SIZE_0) * sizeof(int));

        for(int k = 0; k < n_processes; k++){
            int size = k == n_processes -1 ? NODE_SIZE_0 : NODE_SIZE;
            for (int i = 0; i < size; i++) {
                x[i] = rand() % MAX_VALUE;
               // printf("Gerando x[%d] = %d ", i, x[i]);
            }
           // printf("\n");
            for (int numero_processo = 1; numero_processo < n_processes; numero_processo++) {
                    for (int j = 0; j < NODE_SIZE/n_processes; j++) {
                        aux[j] = x[j * n_processes + numero_processo];
                        if(numero_processo == 1){
                           printf("aux[%d] = %d ", j, aux[j]);
                        }

                    }
                    MPI_Send(aux, NODE_SIZE/n_processes, MPI_INT, numero_processo, messageTag, MPI_COMM_WORLD);           
            }
            printf("\n");
            messageTag++;

            printf("Atribuindo 0\n");
            for (int j = 0; j < NODE_SIZE/n_processes; j++) {
              //  printf("index %d j %d j*n%d x%d\n", index_x0, j, j*n_processes, x[j*n_processes]);
                x0[index_x0++] = x[j * n_processes];
            }

             for (int j = NODE_SIZE; j < size; j++) {
               // printf("index %d j %d j*n%d x%d\n", index_x0, j, j*n_processes, x[j]);
                x0[index_x0++] = x[j];
            }

        }
        
        
        // messageTag = 0;
        for(int k = 0; k < n_processes; k++){
            int size = k == n_processes ? NODE_SIZE_0 : NODE_SIZE;
            for (int i = 0; i < size; i++) {
                y[i] = rand() % MAX_VALUE;
            }
            for (int numero_processo = 1; numero_processo < n_processes; numero_processo++) {
                    for (int j = 0; j < NODE_SIZE/n_processes; j++) {
                        aux[j] = y[j * n_processes + numero_processo];
                    }
            
                    MPI_Send(aux, NODE_SIZE/n_processes, MPI_INT, numero_processo, messageTag, MPI_COMM_WORLD);   
            }
            messageTag++;
            for (int j = 0; j < size/n_processes; j++) {
                y0[index_y0++] = y[j * n_processes];
            }
        }

        // messageTag = 0;
        
        for(int k = 0; k < n_processes; k++){
            int size = k == n_processes ? NODE_SIZE_0 : NODE_SIZE;
            for (int i = 0; i < size; i++) {
                z[i] = rand() % MAX_VALUE;
            }
            for (int numero_processo = 1; numero_processo < n_processes; numero_processo++) {
                for (int j = 0; j < NODE_SIZE/n_processes; j++) {
                    aux[j] = z[j * n_processes + numero_processo];
                }
                MPI_Send(aux, NODE_SIZE/n_processes, MPI_INT, numero_processo, messageTag, MPI_COMM_WORLD);
            }
            messageTag++;
            for (int j = 0; j < size/n_processes; j++) {
                z0[index_z0++] = z[j * n_processes];
            }
        }
  
        
  
        memcpy(x, x0, NODE_SIZE_0 * sizeof *x);
        memcpy(y, y0, NODE_SIZE_0 * sizeof *y);
        memcpy(z, z0, NODE_SIZE_0 * sizeof *z);
     
    } else {
        // TODO: o tamanho nodesize no receive ta certo?
        int mensageTag = 0;
        for (int i = 0; i < n_processes; i++) {
            MPI_Recv(x + i * NODE_SIZE/n_processes, NODE_SIZE, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
        }

        for (int i = 0; i < n_processes; i++) {
            MPI_Recv(y+ i * NODE_SIZE/n_processes, NODE_SIZE, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
        }

        for (int i = 0; i < n_processes; i++) {
            MPI_Recv(z + i * NODE_SIZE/n_processes, NODE_SIZE, MPI_INT, MASTER_PROCESS, mensageTag++, MPI_COMM_WORLD, &status);
        }

    }
    if(myrank == 0){
        printf("RECEBIDO %d\n", myrank);
         for(int i = 0; i < NODE_SIZE_0; i++){
            // printf("x[%d] = %d ", i, x[i]);
         }
         printf("\n");
    }
    // 
     
    //TODO: n deveria ser pares?

    for (int i = 0; i < n_processes; i++) {
        
        int *x_recebido = (int *) malloc(NODE_SIZE_0 * sizeof(int));
        int *y_recebido = (int *) malloc(NODE_SIZE_0 * sizeof(int));
        int *z_recebido = (int *) malloc(NODE_SIZE_0 * sizeof(int));
        int size_i = i == 0 ? NODE_SIZE_0 : NODE_SIZE;

        if(myrank == i){
           // printf("My rank %d faz broadcast\n", i);
            MPI_Bcast(x, size_i, MPI_INT, myrank, MPI_COMM_WORLD);
            MPI_Bcast(y, size_i, MPI_INT, myrank, MPI_COMM_WORLD);
            MPI_Bcast(z, size_i, MPI_INT, myrank, MPI_COMM_WORLD);
            memcpy(x_recebido, x, size_i * sizeof(int));
            memcpy(y_recebido, y, size_i * sizeof(int));
            memcpy(z_recebido, z, size_i * sizeof(int));
        }else{
            MPI_Bcast(x_recebido, size_i, MPI_INT, i, MPI_COMM_WORLD);
            MPI_Bcast(y_recebido, size_i, MPI_INT, i, MPI_COMM_WORLD);
            MPI_Bcast(z_recebido, size_i, MPI_INT, i, MPI_COMM_WORLD);
        }
            int size_j = myrank == 0 ? NODE_SIZE_0 : NODE_SIZE;

            for (int j = 0; j < size_j; j++) { // posição j -> menor distancia em relação ao vetor que tá no processo i
                // int minManhattanLocal_parcial = INT_MAX, maxManhattanLocal_parcial = INT_MIN;
                // double minEuclideanLocal_parcial = INT_MAX, maxEuclideanLocal_parcial = INT_MIN;

            
                for (int k = 0; k < size_i; k++) {
                    if ((j * n_processes + myrank) < (k * n_processes + i) || (k > NODE_SIZE && j < k)) {
                       
                        // K: dados do processo que ta chegando por broadcast
                        // J: dados do processo atual
                        int manhattanDistance = manhattan(x[j], y[j], z[j], x_recebido[k], y_recebido[k], z_recebido[k]);
                      //  if(myrank == 0)
                        //    printf("to no myrank %d e olhando %d calculei a distancia de %d % d = %d\n", myrank, i, j, k, manhattanDistance);
                        if(myrank == 3){
                            // printf("Manhattan %d\n", manhattanDistance);
                            // printf("J = %d | k = %d \n", j, k);
                            // printf("Real j = %d | Real k = %d\n",(j * n_processes + myrank) , (k * n_processes + i));
                            // printf("x %d x_r %d\n", x[j], x_recebido[k]);
                            // printf("y %d y_r %d\n", y[j], y_recebido[k]);
                            // printf("z %d z_r %d\n", z[j], z_recebido[k]);

                        }
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
            } // end 
        free(x_recebido);
        free(y_recebido);
        free(z_recebido);
    } // end i
    int size = myrank == 0 ? NODE_SIZE_0 : NODE_SIZE;

   // printf("Printando distancias max do myrank %d\n", myrank);
    for(int i = 0; i < size; i++){
        //printf("%d ", dist_max_parcial_manhattan[i]);
    }
    printf("\n");
    // printf("Printando distancias min do myrank %d\n", myrank);
    // for(int i = 0; i < NODE_SIZE; i++){
    //     printf("%d ", dist_min_parcial_manhattan[i]);
    // }
    // printf("\n");
    
    int minLocalAbsolutoManhattan = INT_MAX, maxLocalAbsolutoManhattan = INT_MIN, sumMinLocalManhattan = 0, sumMaxLocalManhattan = 0; 
    double minLocalAbsolutoEuclidiana = INT_MAX, maxLocalAbsolutoEuclidiana = INT_MIN, sumMinLocalEuclidiana = 0, sumMaxLocalEuclidiana =0; 
   // int size = myrank == 0 ? NODE_SIZE_0 : NODE_SIZE;
    printf("my rank %d size %d\n", myrank, size);
    #pragma omp parallel for \
        reduction(min: minLocalAbsolutoManhattan, minLocalAbsolutoEuclidiana) \
        reduction(max: maxLocalAbsolutoManhattan, maxLocalAbsolutoEuclidiana) \
        reduction(+: sumMinLocalManhattan, sumMaxLocalManhattan, sumMinLocalEuclidiana, sumMaxLocalEuclidiana)
    for(int i = 0; i < size; i++){
       if (dist_min_parcial_manhattan[i] < minLocalAbsolutoManhattan) {
            minLocalAbsolutoManhattan = dist_min_parcial_manhattan[i];
        }

        if (dist_max_parcial_manhattan[i] > maxLocalAbsolutoManhattan) {
            maxLocalAbsolutoManhattan = dist_max_parcial_manhattan[i];
        }

        if (dist_min_parcial_euclidiana[i] < minLocalAbsolutoEuclidiana) {
            minLocalAbsolutoEuclidiana = dist_min_parcial_euclidiana[i];
        }

        if (dist_max_parcial_euclidiana[i] > maxLocalAbsolutoEuclidiana) {
            maxLocalAbsolutoEuclidiana = dist_max_parcial_euclidiana[i] ;
        }

        // atualiza a soma dos valores máximos e mínimos
        if (dist_min_parcial_manhattan[i] != INT_MAX) {
            sumMaxLocalManhattan += dist_max_parcial_manhattan[i];
            sumMinLocalEuclidiana += dist_min_parcial_euclidiana[i];
            sumMinLocalManhattan += dist_min_parcial_manhattan[i];
            sumMaxLocalEuclidiana += dist_max_parcial_euclidiana[i];
        }
    }
    // fazer reduce entre todos os processos, pois temos 4 vetores de distancias
    // 

    MPI_Reduce(&minLocalAbsolutoManhattan, &minManhattanGlobal, 1, MPI_INT, MPI_MIN, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&maxLocalAbsolutoManhattan, &maxManhattanGlobal, 1, MPI_INT, MPI_MAX, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&minLocalAbsolutoEuclidiana, &minEuclideanGlobal, 1, MPI_DOUBLE, MPI_MIN, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&maxLocalAbsolutoEuclidiana, &maxEuclideanGlobal, 1, MPI_DOUBLE, MPI_MAX, MASTER_PROCESS, MPI_COMM_WORLD);

    MPI_Reduce(&sumMinLocalManhattan, &somaMinManhattan, 1, MPI_INT, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&sumMaxLocalManhattan, &somaMaxManhattan, 1, MPI_INT, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&sumMinLocalEuclidiana, &somaMinEuclidiana, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Reduce(&sumMaxLocalEuclidiana, &somaMaxEuclidiana, 1, MPI_DOUBLE, MPI_SUM, MASTER_PROCESS, MPI_COMM_WORLD);


        //
    if(myrank == 0){
        printf("Distância de Manhattan mínima: %d (soma min: %d) e máxima: %d (soma max: %d).\n", minManhattanGlobal, somaMinManhattan, maxManhattanGlobal, somaMaxManhattan);
        printf("Distância Euclidiana mínima: %.2lf (soma min: %.2lf) e máxima: %.2lf (soma max: %.2lf).\n", minEuclideanGlobal, somaMinEuclidiana, maxEuclideanGlobal, somaMaxEuclidiana);
    }
    free(x);
    free(y);
    free(z);
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