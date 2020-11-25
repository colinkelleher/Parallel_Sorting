// ==================================
// Simple Sort
// ==================================

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"



double * merge_array(int n, double * a, int m, double * b);
void     merge_sort(int n, double * a);
void     swap (double * a, double * b);

int MPI_Array_sort(int n, double * a, int root, MPI_Comm comm);

int main (int argc, char *argv[])
{
   //setup
    int rank, size;
    int n = 10000000, i, j, k, x, q, l, shell, pair, *nr;
    double m = 10.0;
    double * scattered_array;
    double *array, processorTime;

    // Init + rank + size
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if( rank == 0 )
    {

       //initialise the array with random values, then scatter to all processors
       array = (double *) calloc( n, sizeof(double) );
       srand( ((unsigned)time(NULL)+rank) );

       for( i = 0; i < n; i++ )
       {
          array[i]=((double)rand()/RAND_MAX)*m;
         //  printf( "Input: %f \n", array[i] );
       }

    }
   processorTime = MPI_Wtime( );

    MPI_Array_sort(n, array, 0, MPI_COMM_WORLD);

   processorTime = MPI_Wtime( ) - processorTime;
   printf( "Processor %d runs for %lf sec\n", rank, processorTime );

    if( rank == 0 )
    {

      // print array
      // for(int i=0; i<n; i++)printf("Output: %f\n",array[i]);

    }

    

    MPI_Finalize();

}

// MPI Functions

int MPI_Array_sort(int n, double * a, int root, MPI_Comm comm){


    // get rank and size of comm
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

   double commT = 0, compT = 0, time;

    // allocate localA
    double * localA = (double *)calloc(n/size, sizeof(double));

    // scatter a to localA
    int error = MPI_Scatter(a , n/size , MPI_DOUBLE , localA ,n/size , MPI_DOUBLE , root , comm);
    if(error != MPI_SUCCESS) return error;


    // sort localA
    time = MPI_Wtime();
    merge_sort(n/size, localA);
    // bubble_sort(n/size, localA)
    time = MPI_Wtime()-time;
    compT += time;

    // gather localA to a
    time = MPI_Wtime();
    error = MPI_Gather( localA , n/size , MPI_DOUBLE , a , n/size , MPI_DOUBLE , root , comm);
    if(error != MPI_SUCCESS) return error;
    time = MPI_Wtime()-time;
    commT += time;

    // if rank 0 then restore the order
    time = MPI_Wtime();
    if(rank == 0){
        // perform size-1 merges
        for(int i = 1; i<size; i++){
            double * tmp = merge_array(i*n/size, a, n/size, a+i*n/size);
            for(int j=0; j<(i+1)*n/size; j++)a[j] = tmp[j];
        }
    }
    time = MPI_Wtime()-time;
    compT += time;

    printf("Processor %d commT = %lf compT = %lf\n", rank, commT, compT);
    return MPI_SUCCESS;
}



// function to merge the array a with n elements with the array b with m elements
// function returns the nerged array

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

// function to merge sort the array a with n elements

void merge_sort(int n, double * a){

   double * c;
   int i;

   if (n<=1) return;

   if(n==2) {

      if(a[0]>a[1])swap(&a[0],&a[1]);
      return;
   }



   merge_sort(n/2,a);merge_sort(n-n/2,a+n/2);

   c=merge_array(n/2,a,n-n/2,a+n/2);

   for(i=0;i<n;i++)a[i]=c[i];

return;
}

// swap two doubles
void swap (double * a, double * b){

   double temp;

   temp=*a;*a=*b;*b=temp;

}