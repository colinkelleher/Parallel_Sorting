// ==================================
// Bucket Sort
// ==================================

// include headers
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
 
double * merge_array(int n, double * a, int m, double * b);
void merge_sort(int n, double * a);
void swap (double * a, double * b);
 
 
 
// int MPI_Sort_ranking(int n, double * a, double max, int root, MPI_Comm comm);
int  MPI_Sort_bucket(int n, double * a, double max, int root, MPI_Comm comm);

// ------------------------------
 
int main (int argc, char *argv[]){
 
	int rank, size;
 
	int n = 1000000, q, l, i, j, k, x, *nr;
	double m = 10.0;
	double *a, *b;
 
	MPI_Status status;
 
	MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
 
	a = (double *) calloc(n,sizeof(double));
	b = (double *) calloc(n,sizeof(double));
 
	if( rank == 0 )
	{
 
	   //initialise the array with random values, then scatter to all processors
	   srand( ((unsigned)time(NULL)+rank) );
 
	   for( i = 0; i < n; i++ )
	   {
	      a[i]=((double)rand()/RAND_MAX)*m;
	      //printf( "Initial: %f\n", a[i] );
	   }
 
	}
 
   double time = MPI_Wtime();
	MPI_Sort_bucket(n, a, m, 0, MPI_COMM_WORLD);
   time = MPI_Wtime()-time;
 
    printf("Processor %d runs for %lf\n", rank, time);
 
	if( rank == 0 )
	{
	   for( i = 0; i < n; i++ )
	   {
	      //printf( "Output : %f\n", a[i] );
	   }
	}
	MPI_Finalize();
 
}
 
int MPI_Sort_bucket(int n, double * a, double max, int root, MPI_Comm comm)
{
 
    // find rank and size
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
 
    double commT = 0, compT = 0, time;
 
    // allocate the extra memory / arrays needed
    double * bucket = (double *)calloc(n, sizeof(double));
 
    // Brodcast the array to the processor
 
    time = MPI_Wtime();
    int error = MPI_Bcast(a, n, MPI_DOUBLE, root, comm);
    if(error != MPI_SUCCESS)return error;
    time = MPI_Wtime()-time;
    commT += time;
 
    // Scan the array to find bucket
 
    time = MPI_Wtime();
    int count = 0;
    for(int i=0;i<n;i++){
        if(rank*max/size<=a[i] && a[i]<(rank+1)*max/size)
            bucket[count++] = a[i];
    }
 
    // Sort bucket
    merge_sort(count, bucket);
    time = MPI_Wtime()-time;
    compT += time;
 
    // Gatherv bucket
    int * counts = (int *)calloc(size, sizeof(int));
    int * displs = (int *)calloc(size, sizeof(int));
 
    time = MPI_Wtime();
    MPI_Gather(&count, 1, MPI_INT, counts, 1, MPI_INT, root, comm);
    time = MPI_Wtime()-time;
    commT += time;
 
    time = MPI_Wtime();
    displs[0] = 0;
    for(int i=1;i<size;i++)
        displs[i] = displs[i-1] + counts[i-1];
    time = MPI_Wtime()-time;
    compT += time;
 
    time = MPI_Wtime();
    error = MPI_Gatherv(bucket, count, MPI_DOUBLE, a, counts, displs, MPI_DOUBLE, root, comm);
    if(error != MPI_SUCCESS)return error;
    time = MPI_Wtime()-time;
    commT += time;
 
    printf("Processor %d commT = %lf compT = %lf\n", rank, commT, compT);
 
    return MPI_SUCCESS;
 
}
 
 
void  MPI_Sort_ranking(int n, double * a, double max, int root, MPI_Comm comm)
{
 
	// find rank and size
 
	// allocate the extra memory / arrays needed
 
	// Brodcast the array to the processor
 
	// P rank generates an array ranking with ranking[i] is the rank of a[i+rank*n/size] in the array
 
	// Gather the array ranking to finalRanking
 
	// if processor 0 then restore the order in the array b and move b back to a
 
 
}
 
 
// ------------------------------------------------------------
//
// these functions deal with sorting and merging
//
 
double * merge_array(int n, double * a, int m, double * b){
   int i,j,k;
   double * c = (double *) calloc(n+m, sizeof(double));
 
   for(i=j=k=0;(i<n)&&(j<m);)
      if(a[i]<=b[j])c[k++]=a[i++];
      else c[k++]=b[j++];
      if(i==n)for(;j<m;)c[k++]=b[j++];
      else for(;i<n;)c[k++]=a[i++];
   return c;
}
 
void merge_sort(int n, double * a){
   double * c;
   int i;
   if (n<=1) return;
   if(n==2) {
      if(a[0]>a[1])swap(&a[0],&a[1]);
      return;
   }
   merge_sort(n/2,a);
   merge_sort(n-n/2,a+n/2);
   c=merge_array(n/2,a,n-n/2,a+n/2);
   for(i=0;i<n;i++)a[i]=c[i];
}
 
void swap (double * a, double * b){
 
   double temp;
 

 
}