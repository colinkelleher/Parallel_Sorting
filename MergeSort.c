// ==================================
// Merge Sort
// ==================================
#include <time.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void merge_sort(int, double *);
void merge_array(int, double *, int, double *);
void swap (double *, double *);

int isActive(int rank, int p, int l);
int isSender(int rank, int p, int l);
int isReciever(int rank, int p, int l);

int main (int argc, char *argv[]) {
	int rank, size;
    int n = 1000000;
	int q, l, i, j, k, x, *nr;
	double *a, *b, m = 10.0, processorTime;

    MPI_Status status;
    MPI_Init(&argc, &argv);

   	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   	MPI_Comm_size(MPI_COMM_WORLD, &size);

	a = (double *) calloc(n,sizeof(double));
    b = (double *) calloc(n,sizeof(double));

    processorTime = MPI_Wtime( );

	if(rank == 0) {
		srand(((unsigned)time(NULL)+rank)); 		// Initialise the array with random values
		for(x = 0; x < n; x++){
			a[x]=((double)rand()/RAND_MAX)*m;
		}
	}

	q = 0;
	while(size > pow(2,q)){
        q++;
    }

    double commT = 0, compT = 0, time;
    
    // Stage 1	
	for(l = 0; l <= q; l++){
		if( isReciever(rank, size, l) && l > 0){
            time = MPI_Wtime();
			MPI_Recv(a, n/(int)pow(2, l), MPI_DOUBLE, rank - size/(int)pow(2, l), 0, MPI_COMM_WORLD, NULL);
            time = MPI_Wtime()-time;
            commT += time;
        }
		if( isActive(rank, size, l) && l < q ){
			MPI_Send(&a[n/(int)pow(2, l+1)], n/(int)pow(2, l+1), MPI_DOUBLE, rank+size/(int)pow(2, l+1), 0, MPI_COMM_WORLD);
        }
	}

    // Stage 2
    time = MPI_Wtime();
	merge_sort(n/size, a);
    time = MPI_Wtime()-time;
    compT += time;
	
    // Stage 3
	for(l = q; l >= 0; l--){
		if( isActive(rank, size, l) && l < q ){
            time = MPI_Wtime();
			x = n/(int)pow(2, l+1);
			MPI_Recv(&a[x], x, MPI_DOUBLE, rank+size/(int)pow(2, l+1), 0, MPI_COMM_WORLD, NULL); 
			    time = MPI_Wtime()-time;
                commT += time;



            time = MPI_Wtime();
			merge_array(x, a, x, &a[x]);
            time = MPI_Wtime()-time;
            compT += time;
		}
		if( isSender(rank, size, l) && l > 0 )
			MPI_Send(a, n/(int)pow(2, l), MPI_DOUBLE, rank-size/(int)pow(2, l), 0, MPI_COMM_WORLD);
	}
	processorTime = MPI_Wtime( ) - processorTime;
    printf( "Processor %d runs for %lf sec\n", rank, processorTime );
    printf("Processor %d commT = %lf compT = %lf\n", rank, commT, compT);
	MPI_Finalize();

    return 0;
}

int isActive(int rank, int p, int l) {
	if((rank % (int)(p/pow(2,l))) == 0 ) 
		return 1;
	return 0;
}

int isSender(int rank, int p, int l) {
   return isReciever( rank, p , l );
}

int isReciever(int rank, int p, int l) {
	if(!isActive(rank, p, l))  return 0;
	int result = rank / ((int)(p/pow(2,l)) );
	if( result % 2 == 1) return 1;
	return 0;
}


// as before

void merge_array(int n, double * a, int m, double * b){
	int i,j,k;
	double * c = (double *) calloc(n+m, sizeof(double));
   
	for(i=j=k=0; (i<n)&&(j<m);)
		if(a[i] <= b[j]) c[k++] = a[i++];
		else c[k++] = b[j++];
	if(i==n) 
		for(;j<m;) c[k++] = b[j++];
	else
		for(;i<n;) c[k++] = a[i++];
	for(i=0; i<n+m; i++)
		a[i] = c[i];
}

void merge_sort(int n, double * a){
	double * c;
	int i;
	if(n<=1) return;
	if(n==2) {
		if(a[0] > a[1]) 
			swap(&a[0],&a[1]);
      return;
	}
	merge_sort(n/2, a);
	merge_sort(n-n/2, a+n/2);
	merge_array(n/2, a, n-n/2, a+n/2);
}

void swap (double * a, double * b){
   double temp;
   temp=*a;*a=*b;*b=temp;
}
