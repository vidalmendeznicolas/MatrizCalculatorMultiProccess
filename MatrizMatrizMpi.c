#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include <math.h>

// mpicc MatrizMatrizMpi.c -o MatrizMatrizMpi -lm
// mpirun -np 4 ./MatrizMatrizMpi 4 6 4 2 1
// mpirun -np 4 ./MatrizMatrizMpi 4 6 4 1 1
// mpirun --oversubscribe -np 4 ./MatrizMatrizMpi 2 2 4 1 1 
//        procs     ->  m n k alfa test
// 0 ===> Creará comunicadores con split
// 1 ===> Creará subtopologias virtuales
#define CART 0 
#define LIM 600

/* matrizMatrizLocal multiplica una matriz por un vector en memoria local */
void matrizMatrizLocal(float alpha, float *A, int filasA, int colAfilB, int colB, float *B, float *y, int rank) {
	// int o=0;
	// printf("Entro a multiplicar con proceso %d \n", rank);
	// 	for(o=0;o<filasA*colAfilB;o++) {
	// 		printf("A es %f con rank %d \n", A[o], rank);
	// 	}
	// 	for(o=0;o<colAfilB*colB;o++) {
	// 		printf("B es %f con rank %d \n", B[o], rank);
	// 	}
	int i,j,z;	
	for (i=0; i<filasA; i++) {					//filas A -> filas y
      	for (j=0; j<colB; j++) {				//columnas B -> columnas y
			for (z=0; z<colAfilB; z++) { 		//columnas A, filas B 
				//printf("y[%d] += %d*A[%d]*B[%d] // %d += %d*%d*%d rank: %d\n",i*colB+j,(int)alpha,i*colAfilB+z,z*colB+j,(int)y[i*colB+j],(int)alpha,(int)A[i*colAfilB+z],(int)B[z*colB+j],rank);
        		y[i*colB+j] += alpha*A[i*colAfilB+z]*B[z*colB+j];
				//printf("y[%d] += %d*A[%d]*B[%d] // %d += %d*%d*%d rank: %d\n",i*colB+j,(int)alpha,i*colAfilB+z,z*colB+j,(int)y[i*colB+j],(int)alpha,(int)A[i*colAfilB+z],(int)B[z*colB+j],rank);
			}
		}
	}
	return;
}

int main(int argc, char *argv[]) {
    int numprocesos, rank, i, j, m, n, k, test, sqrtnumprocesos, anchoA, altoA, anchoB, altoB, anchoC, altoC, dim;
    int position = 0;
    float alpha;
	char buffer[LIM];

    MPI_Init(&argc, &argv);
    // Determinar el rango del proceso invocado
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Determinar el numero de procesos
    MPI_Comm_size(MPI_COMM_WORLD, &numprocesos);

    dim = sqrt(numprocesos);

    // Proceso 0 lee parámetros de entrada
    // Parámetro 1 -> m - Número de filas de la matriz A
    // Parámetro 2 -> n - Número de Columnas de la matriz B
    // Parámetro 3 -> k - Número de filas de la matriz B y columnas de la matriz A
	// Parámetro 4 -> alpha
    // Parámetro 5 -> booleano que nos indica si se desea imprimir matrices y vectores de entrada y salida
    if(!rank){
        // Seguir solo si es una malla cuadrada
        if ((sqrt(numprocesos) - (int)sqrt(numprocesos))) {
			printf("NUMERO DE PROCESOS INCORRECTO: No se puede generar una malla cuadrada \n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			MPI_Finalize();
			return 0;
		}
        if(argc>5){
	    	m = atoi(argv[1]);
            n = atoi(argv[2]);
	    	k = atoi(argv[3]);
            alpha = atof(argv[4]);
            test = atoi(argv[5]);
			if (m % (int)dim || n % (int)dim || k % (int)dim ) {
				printf("PARAMETROS INCORRECTOS: m,n y k tienen que ser multiplos de sqrt(numprocesos)\n");
				MPI_Abort(MPI_COMM_WORLD, 1);
				MPI_Finalize();
				return 0;
			}
        }
        else{
            printf("NUMERO DE PARAMETROS INCORRECTO\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
            return 0;
        }

    }
	
    //Se empaquetes m,n,k y alpha en el proceso 0
	if (!rank) {
		MPI_Pack(&n, 1, MPI_INT, buffer, LIM, &position, MPI_COMM_WORLD);
		MPI_Pack(&k, 1, MPI_INT, buffer, LIM, &position, MPI_COMM_WORLD);
		MPI_Pack(&m, 1, MPI_INT, buffer, LIM, &position, MPI_COMM_WORLD);
		MPI_Pack(&alpha, 1, MPI_FLOAT, buffer, LIM, &position, MPI_COMM_WORLD);
	}
    // El proceso 0 envía alpha, m, n, k al resto de procesos de forma empaquetada
	MPI_Bcast(&buffer, LIM, MPI_PACKED, 0, MPI_COMM_WORLD);
    //El resto de procesos desempaquetan m, n, k y alpha
	if (rank) {
		MPI_Unpack(buffer, LIM, &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
		MPI_Unpack(buffer, LIM, &position, &k, 1, MPI_INT, MPI_COMM_WORLD);
		MPI_Unpack(buffer, LIM, &position, &m, 1, MPI_INT, MPI_COMM_WORLD);
		MPI_Unpack(buffer, LIM, &position, &alpha, 1, MPI_FLOAT, MPI_COMM_WORLD);
	}

	sqrtnumprocesos = (int)dim;
	anchoA = k/sqrtnumprocesos;
	altoA = m/sqrtnumprocesos;
	altoB = anchoA;
	anchoB = n/sqrtnumprocesos;
	altoC = altoA;
	anchoC = anchoB;

	MPI_Datatype submatrizA, submatrizB, submatrizC;
	MPI_Comm comFilas, comColumnas;

// definimos los tamaños de los cachitos de la matriz
// el stride es el numero respecto al 0 de la submatriz, que nos indica cual es la siguiente posicion.
    MPI_Type_vector(altoA, anchoA, k, MPI_FLOAT, &submatrizA);
	MPI_Type_vector(altoB, anchoB, n, MPI_FLOAT, &submatrizB);
	MPI_Type_vector(altoC, anchoC, n, MPI_FLOAT, &submatrizC);
	MPI_Type_commit(&submatrizA);
	MPI_Type_commit(&submatrizB);
	MPI_Type_commit(&submatrizC);

    if (CART) {
        int aux1[] = {sqrtnumprocesos,sqrtnumprocesos};
        int aux2[] = {0,0};
        MPI_Comm newComunicador;
        MPI_Cart_create(MPI_COMM_WORLD, 2, aux1, aux2, 0, &newComunicador);
        aux1[0] = 0; aux1[1] = 1;
        MPI_Cart_sub(newComunicador, aux1, &comFilas);
        aux1[0] = 1; aux1[1] = 0;
        MPI_Cart_sub(newComunicador, aux1, &comColumnas);
    }
    else {
        MPI_Comm_split(MPI_COMM_WORLD, rank/sqrtnumprocesos, rank % sqrtnumprocesos, &comFilas);
        MPI_Comm_split(MPI_COMM_WORLD, rank % sqrtnumprocesos, rank/sqrtnumprocesos, &comColumnas);
    }

    float *B;
    float *A;
    float *y;
	float *localA = (float *) malloc(anchoA*altoA*sizeof(float));
	float *localB = (float *) malloc(anchoB*altoB*sizeof(float));

     if(!rank){ // Proceso 0 inicializa la matriz y el vector
        A = (float *) malloc(m*k*sizeof(float));
		B = (float *) malloc(k*n*sizeof(float));
		y = (float *) malloc(m*n*sizeof(float));

        for(i=0; i<m; i++){		//Matriz A
            for(j=0; j<k; j++){
                A[i*k+j] = 1+i+j;
            }
        }

        for(i=0; i<k; i++){		//Matriz B
            for(j=0; j<n; j++){
                B[i*n+j] = 1+i+j;
            }
        }

        if(test){
            printf("\nMatriz A es...\n");
            for(i=0; i<m; i++){
                for(j=0; j<k; j++){
                    printf("%f ", A[i*k+j]);
                }
                printf("\n");
            }
            printf("\nMatriz B es...\n");
            for(i=0; i<k; i++){
                for(j=0; j<n; j++){
                    printf("%f ", B[i*n+j]);
                }
                printf("\n");
            }
		}
    }

   
	double time;
    MPI_Barrier(MPI_COMM_WORLD); // Barrera para garantizar una correcta medida de tiempo
    time = MPI_Wtime();
    // Vamos a dividir las matrices en submatrices y distribuirlas en los procesos
	if (!rank) {
		int indiceA=0; // posición actual en A
        int indiceB=0; // posición actual en B
        // Determinamos el tamaño de cada submatriz y cómo avanzar a la siguiente submatriz en la matriz original.
		int nextbloqueA = anchoA;                    //avanzamos un bloque en la misma fila
		int nextbloqueA2 = k*m/sqrtnumprocesos - nextbloqueA*(sqrtnumprocesos -1);
		int nextbloqueB = anchoB;
		int nextbloqueB2 = k*n/sqrtnumprocesos - nextbloqueB*(sqrtnumprocesos -1);

		for(i=0; i<numprocesos; i++) {
			if (i != 0 && (i % sqrtnumprocesos == 0)) {	//vamos a la sig fila
				indiceA += nextbloqueA2;
				indiceB += nextbloqueB2; 
			} else if (i != 0) {	                    //vamos al sig bloque de la misma fila
				indiceA += nextbloqueA; 
				indiceB += nextbloqueB;
			}
			MPI_Send(&A[indiceA],1, submatrizA, i, i, MPI_COMM_WORLD);
			MPI_Send(&B[indiceB],1, submatrizB, i, i, MPI_COMM_WORLD);
		}
	} 

	MPI_Recv(localA, anchoA*altoA, MPI_FLOAT, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(localB, anchoB*altoB, MPI_FLOAT, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /****IMPLEMENTACIÓN DEL ALGORITMO SUMMA******/	
	// Se crean tres buferes para almacenar las submatrices locales y el resultado de la multiplicación
	float *bufA = (float *) malloc(anchoA*altoA*sizeof(float));
	float *bufB = (float *) malloc(anchoB*altoB*sizeof(float));
	float *localC = (float *) malloc(anchoC*altoC*sizeof(float));

    //Todos los procesos inicia su submatriz C a 0
	for(i=0;i<anchoC*altoC;i++) {
		localC[i]=0.0;
	}

	for(i=0;i<sqrtnumprocesos;i++) {
        // Se copian los valores de las submatrices locales en los bufers para ser enviados
		for(j=0;j<altoA*anchoA;j++) {
			//printf("localA es %f con rank %d \n", localA[j], rank);
			bufA[j] = localA[j];
		}
		for(j=0;j<altoB*anchoB;j++) {
			//printf("localB es %f con rank %d \n", localB[j], rank);
			bufB[j] = localB[j];
		}
        // Se transmiten los buferes a los procesos porque cada proceso necesita una copia de la submatriz completa para
        // su parte de la multiplicación de matrices
		MPI_Bcast(bufA, anchoA*altoA, MPI_FLOAT, i, comFilas);
		MPI_Bcast(bufB, anchoB*altoB, MPI_FLOAT, i, comColumnas);
        //Se multiplica la submatriz por la matriz y almacena el resultado en localC
		//printf("iteracion %d del proc %d \n", i, rank);
		matrizMatrizLocal(alpha, bufA, altoA, anchoA, anchoB, bufB, localC,rank); 
	}
    // Enviamos localC al proceso 0
	MPI_Send(localC, anchoC*altoC, MPI_FLOAT, 0, rank, MPI_COMM_WORLD);
    time = MPI_Wtime() - time;
    if(!rank){
		int proximobloque = anchoC;                     //avanzamos un bloque en la misma fila
		int proximobloque2 = m*n/sqrtnumprocesos - proximobloque*(sqrtnumprocesos -1);
		int indice=0;
		for(i=0;i<numprocesos;i++) {
			if ((i % sqrtnumprocesos == 0) && i != 0) {	//vamos a la sig fila
				indice += proximobloque2;
			} else if (i != 0) {	                    //vamos al sig bloque de la misma fila
				indice += proximobloque; 
			}
			MPI_Recv(&y[indice],1, submatrizC, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	
        if(test){
            printf("\nAl final la matriz resultado \"y\" es...\n");
            for(i=0; i<m; i++){
				for(j=0; j<n; j++) {
                	printf("%f ", y[i*n +j]);
				}
				printf("\n");
            }
            printf("\n");

            // Solo el proceso 0 calcula el producto y compara los resultados del programa secuencial con el paralelo
            float *test = (float *) malloc(m*n*sizeof(float));
            for(i=0; i<m; i++){
				for(j=0; j<n; j++) {
					test[i*n +j] = 0.0;
				}
            }

			matrizMatrizLocal(alpha, A, m, k, n, B, test,0); 
			
            int err = 0;
            for(i=0; i<m; i++){
				for(j=0; j<n; j++) {
					if(test[i*n+j] != y[i*n+j]){
						err++;
						printf("\n Error en la posicion %d porque %f != %f", i*n+j, y[i*n+j], test[i*n+j]);
					}
				}	
			}
            printf("\n%d errores en el producto matriz A %dx%d matriz B %dx%d con %d procesos\n", err, m,k,k,n,numprocesos);
            free(test);
        }

        free(A);
		free(B);
        free(y);
    }

    free(localA);
	free(localB);
    free(localC);
	free(bufA);
	free(bufB);

    // Barrera para que no se mezcle la impresión del tiempo con la de los resultados
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("\nProducto matriz vector con dimension %d y %d procesos: Tiempo de ejecucion del proceso %d fue %lf\n", n, numprocesos, rank, time);

    MPI_Finalize();
    return 0;
}
