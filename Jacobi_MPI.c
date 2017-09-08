#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include<time.h>
#include<mpi.h>
#include<math.h>

int main(int argc, char *argv[]){

		int size,rank;
							
		MPI_Init(&argc,&argv);
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		MPI_Status status;
		double inizio = MPI_Wtime();
		
		srand(7);
							
		int n = 12000;
		int tag =  0;
		int iter = 0;
////////////////////////////////////////////////////////////////////////////master
		if(rank == 0){
		
			double *matrix;
			matrix=(double *)malloc(sizeof(double)* n * n);

			int i,j;
			double diffnorm_finale = 0;
			int begin_row_slave;
		
			for(i = 0;i < n;i++){//riempimento della matrice
				for(j=0;j < n;j++){
					if(i == 0 || i == n-1)
						matrix[(i*n)+j] = -1;
					else
						matrix[(i*n)+j] = rand() % 10;
											
				}
			}	
			
	/*			for(i = 0; i < n ; i++){//stampo matrice
					for(j = 0;j < n;j++){					
						printf("%.2f ",matrix[(i*n)+j]);
					}
						printf("\n");
				}
				printf("\n");*/
////////////////////////////////////////////////size = 1				
			if(size == 1){
			
				double *matrix_support;
				matrix_support=(double *)malloc(sizeof(double) * n * n);
				
				while(iter < 100){
							
					for(i = 0;i < n;i++)//riempimento della matrice supporto
						for(j=0;j < n;j++){
							if(i == 0 || i == (n-1) || j == 0 || j == (n-1))
								matrix_support[(i*n)+j] = matrix[(i*n)+j];
							else 
								matrix_support[(i * n) + j] = (matrix[(i * n) + j - 1] + matrix[(i * n) + j + 1]+ matrix[((i + 1)*n) + j] + matrix[((i - 1) * n) + j])/4;
							diffnorm_finale = diffnorm_finale + (matrix_support[(i*n)+j] - matrix[(i*n)+j]) * (matrix_support[(i*n)+j] - matrix[(i*n)+j]);
					}
					
				for(i = 0;i < n;i++){//riempimento della matrice
					for(j=0;j < n;j++){
						matrix[(i*n)+j] = matrix_support[(i*n)+j] ;
				}
			}
					
					iter++;
					diffnorm_finale = sqrt(diffnorm_finale);
				}
				
			
			}//fine
///////////////////////////////////////////////////////////////////////////////size > 2
		else{
		
			int temp;
			int resto,rest2;
			int begin,tot_rows,master_row;
			int nrighe = n/size;
			double diffnorm_slave =0;
		
				temp = 0;
			
			
				resto = n % size;
				
				if(resto > 0)
					rest2 = resto-1;
			
				for(i = 0;i < size;i++){//mando righe
								
					if(resto > 0){
					
						if(i == 0){				
							tot_rows = (nrighe+1) +1;
							
							temp =temp + tot_rows -1;
							begin=0;
							resto--;
						}
						else{
							begin = temp-1;
							tot_rows = (nrighe + 2) + 1;
						
							temp = temp + tot_rows - 2;
							resto--;							
						}
					}
					else{
						if(i == 0){
							begin = 0;
							tot_rows = nrighe+1;
							temp = temp + tot_rows -1;
						}
						else if(i == (size-1)){
							begin = temp-1;
							tot_rows = nrighe+1;
							temp = temp + tot_rows - 2;
						}
						else{
							begin = temp-1;
							tot_rows = nrighe + 2;
							temp = temp + tot_rows - 2;
						}
					}
					if(i == 0)
						master_row = tot_rows;
					else{
					MPI_Send(&tot_rows,1,MPI_INT,i,tag,MPI_COMM_WORLD); // mando num righe
					MPI_Send(&matrix[begin * n],n * tot_rows,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);}
								
				}//fine mando righe
				
				double *matrix_support;
				matrix_support=(double *)malloc(sizeof(double) * master_row * n);
				int diffnorm_master=0;
				//printf("\ninzia il lavoro del master\n");
				while(iter < 100){
					
					/*for(i = 0; i < master_row ; i++){//stampo matrice
					for(j = 0;j < n;j++){					
						printf("%.2f ",matrix[(i*n)+j]);
					}
						printf("\n");
				}
				printf("\n");*/
					
					for(i = 0;i < master_row;i++)//riempimento della matrice supporto
						for(j=0;j < n;j++){
							if(i == 0 || i == (master_row-1) || j == 0 || j == (n-1))
								matrix_support[(i*n)+j] = matrix[(i*n)+j];
							else 
								matrix_support[(i * n) + j] = (matrix[(i * n) + j - 1] + matrix[(i * n) + j + 1]+ matrix[((i + 1)*n) + j] + matrix[((i - 1) * n) + j])/4;
							diffnorm_master = diffnorm_master + (matrix_support[(i*n)+j] - matrix[(i*n)+j]) * (matrix_support[(i*n)+j] - matrix[(i*n)+j]);
					}
					diffnorm_master = sqrt(diffnorm_master);
					for(i = 0;i < master_row;i++)//riempimento della matrice
						for(j=0;j < n;j++)
							matrix[(i*n)+j] = matrix_support[(i*n)+j] ;
							
					//peppe
					
					//send
					//	printf("\nultima riga mandata:\n");	
					////	for(j = 0;j < n;j++)					
						//	printf("%.2f ",matrix[((master_row-2)*n)+j]);
						////printf("\n");
							MPI_Send(&matrix[(master_row-2)*n],n,MPI_DOUBLE,rank+1,tag,MPI_COMM_WORLD);
					//		printf("\n");
							
						//recv
				//	printf("\nultima riga ricevuta:\n");
						MPI_Recv(&matrix[(master_row-1)*n],n,MPI_DOUBLE,rank+1,tag,MPI_COMM_WORLD,&status);
					//	for(j = 0;j < n;j++)					
					//	printf("%.2f ",matrix[((master_row-1)*n)+j]);
					//	printf("\n\n");
						
					iter++;
					
				}
				diffnorm_finale =diffnorm_finale + diffnorm_master;
				
				for(i = 1;i < size;i++){
					if(i == 1)
						begin_row_slave = (master_row-1)*n;
						//printf("rest2 = %d\n",rest2);
					if(rest2 > 0){
					///	printf("\nvaloew :%d * %d\n",begin_row_slave,(nrighe+1));
						MPI_Recv(&matrix[begin_row_slave],(nrighe+1)*n,MPI_DOUBLE,i,tag,MPI_COMM_WORLD,&status);
						
						begin_row_slave =begin_row_slave + (nrighe+1)*n;
						
						rest2--;
					}
					else{
						//printf("%d * %d\n",begin_row_slave,n);
						MPI_Recv(&matrix[begin_row_slave],(nrighe)*n,MPI_DOUBLE,i,tag,MPI_COMM_WORLD,&status);
						begin_row_slave =begin_row_slave +(nrighe)*n;
					}

								
					MPI_Recv(&diffnorm_slave,1,MPI_DOUBLE,i,tag,MPI_COMM_WORLD,&status);
					diffnorm_finale = diffnorm_finale + diffnorm_slave;
				}
				
				/*printf("\nmatrice finale\n");
				for(i = 0; i < n ; i++){//stampo matrice
					for(j = 0;j < n;j++){					
						printf("%.2f ",matrix[(i*n)+j]);
						fflush(stdout);
					}
						printf("\n");
						fflush(stdout);
				}*/
								
				//diffnorm_finale = sqrt(diffnorm_finale);
				//printf("diffnorm finale = %e\n",diffnorm_finale);
				iter++;
			
		
		}				
	
			double fine = MPI_Wtime();
	printf("tempo : %f\n",fine-inizio);				
		}//fine master
//////////////////////////////////////////////////////////////////////inizio slave
		else{
			int righe;
			double diffnorm_slave;
			int i,j;
			MPI_Recv(&righe,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
			//printf("righe = %d\n",righe);
				
			double *matrix_slave;
			matrix_slave=(double *)malloc(sizeof(double) * righe * n);
								
			double *matrix_output;
			matrix_output=(double *)malloc(sizeof(double) * righe * n);

			MPI_Recv(matrix_slave,n*righe,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);
			
				int first_row = 1; //prenditi la prima riga
				int last_row  = righe-1; //prenditi l'ultima riga
		
				if(rank == (size-1))
					last_row = righe;
					
			while(iter < 100){
			
			//	printf("iter = %d\n",iter);
			//	fflush(stdout);
			/*	for(i = 0; i < righe ; i++){//stampo matrice
					for(j = 0;j < n;j++){					
						printf("%.2f ",matrix_slave[(i*n)+j]);
					}
						printf("\n");
				}
				printf("\n");*/
				int i,j,k;
				for(i = 0;i < righe;i++){//riempimento della matrice output
					for(j=0;j < n;j++){
						if(i == 0 || i == (righe-1) || j == 0 || j == (n-1))
							matrix_output[(i*n)+j] = matrix_slave[(i*n)+j];
						else 
							matrix_output[(i * n) + j] = (matrix_slave[(i * n) + j - 1] + matrix_slave[(i * n) + j + 1]+ matrix_slave[((i + 1)*n) + j] + matrix_slave[((i - 1) * n) + j])/4;
						diffnorm_slave = diffnorm_slave + (matrix_output[(i*n)+j] - matrix_slave[(i*n)+j]) * (matrix_output[(i*n)+j] - matrix_slave[(i*n)+j]);
					}
				}
				
			/*	for(i = 0; i < righe ; i++){//stampo matrice
					for(j = 0;j < n;j++){					
						printf("%.2f ",matrix_output[(i*n)+j]);
					}
						printf("\n");
				}
				printf("\n");*/
				
				diffnorm_slave = sqrt(diffnorm_slave);
				//printf("diffnorm = %e\n",diffnorm_slave);
				for(i = 0;i < righe;i++){//aggiornamento matrice slave
					for(j=0;j < n;j++){
						matrix_slave[(i*n)+j] = matrix_output[(i*n)+j];	
				}
			}
										
				if(rank != (size-1)){//send
					/*printf("ultima riga mandata:\n");	
					for(j = 0;j < n;j++)					
						printf("%.2f ",matrix_output[((last_row-1)*n)+j]);
					printf("\n");*/
						MPI_Send(&matrix_output[(last_row-1)*n],n,MPI_DOUBLE,rank+1,tag,MPI_COMM_WORLD);}
						
						
					//recv
				//	printf("prima riga ricevuta:\n");
						MPI_Recv(&matrix_slave[0],n,MPI_DOUBLE,rank-1,tag,MPI_COMM_WORLD,&status);	
						/*for(j = 0;j < n;j++)					
							printf("%.2f ",matrix_slave[0+j]);*/
						
						
						
					//send
				//	printf("prima riga mandata:\n");	
				/*	for(j = 0;j < n;j++)					
						printf("%.2f ",matrix_output[n+j]);*/
						MPI_Send(&matrix_slave[n],n,MPI_DOUBLE,rank-1,tag,MPI_COMM_WORLD);
						
						//printf("\n");

					
					//	printf("\n");
				if(rank != (size-1)){//recv
					//printf("ultima riga ricevuta:\n");
						MPI_Recv(&matrix_slave[(last_row)*n],n,MPI_DOUBLE,rank+1,tag,MPI_COMM_WORLD,&status);
					/*	for(j = 0;j < n;j++)					
						printf("%.2f ",matrix_slave[((last_row)*n)+j]);*/
						}
					//	printf("\n");
				//fflush(stdout);
				iter++;
			}
				MPI_Send(&matrix_output[first_row*n],n*(last_row - first_row),MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
				MPI_Send(&diffnorm_slave,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
		}

	MPI_Finalize();
	return 0;
}
