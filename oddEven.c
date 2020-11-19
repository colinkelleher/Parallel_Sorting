// include headers
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <mpi.h>
# include <time.h>

//MPI methods
int MPI_Exchange( int n, double * array, int rank1, int rank2, MPI_Comm comm );
int MPI_Sort_oddeven( int n, double * array, int root, MPI_Comm comm );
int MPI_Is_sorted(int n, double * array, int * answer, int root, MPI_Comm);

//all in previous labs
double * merge( int n, double * array, int m, double * b );
void merge_sort( int n, double * array );
void swap ( double * array, double * b );

// function definitions
int main( int argc, char ** argv ) {

  //setup
  int size, rank, result, i, *answer;
  int n = 10000000;
  double m = 10.0;
  double *array, processorTime;
  MPI_Status status;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  //allocate space for an array of doubles, size n
  array = ( double * ) calloc( n, sizeof( double ) );

  //fills array with random values on root proc
  if( rank == 0 )
  {
    //get random values for array & output for testing
    srand( ( ( unsigned ) time( NULL ) + rank ) );
    for( i = 0; i < n; i++ )
    {
      array[ i ] = ( ( double ) rand( ) / RAND_MAX ) * m;
      //printf( "Initial: %f\n", array[ i ] );
    }
  }
  //get start time for each processor
  processorTime = MPI_Wtime( );

  //MPI_Sort does all the heavy work
  result = MPI_Sort_oddeven( n, array, 0, MPI_COMM_WORLD );
  if( result != MPI_SUCCESS )
  {
    return result;
  }
  //get end time for each processor
  processorTime = MPI_Wtime( ) - processorTime;
  printf( "Processor %d runs for %lf sec\n", rank, processorTime );
  MPI_Finalize( );
}

int MPI_Sort_oddeven( int n, double * array, int root, MPI_Comm comm ) {
        // get rank and size of comm
        int rank, size;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
        double commT = 0, compT = 0, time;
        //allocate space for numElements/numProcessors amount of doubles
        double * local_a = (double *)calloc(n/size, sizeof(double));

        //scatter a to local_a
        time = MPI_Wtime();
        MPI_Scatter(array, n/size, MPI_DOUBLE, local_a, n/size, MPI_DOUBLE, root, comm);
        time = MPI_Wtime() - time;
        commT += time;

        //sort local_a using mergeSort
        time = MPI_Wtime();
        merge_sort(n/size, local_a);
        time = MPI_Wtime()-time;
        compT += time;
        
        //odd-even iterations
        for(int step=0;step < size;step++){
            if((step+rank)%2 == 0){
                if(rank<size-1){
                  time = MPI_Wtime();
                  MPI_Exchange(n/size, local_a, rank, rank+1, comm);
                  time = MPI_Wtime() - time;
                  commT += time;
                }

            }else{
                if(rank>0){
                time = MPI_Wtime();
                MPI_Exchange(n/size, local_a, rank-1, rank, comm);
                time = MPI_Wtime() - time;
                commT += time;}
            }

            MPI_Barrier(comm);
            int answer;
            MPI_Is_sorted(n/size, local_a, &answer, root, comm);
            if(answer == 1) break;
        }

        //gather local_a
        time = MPI_Wtime();
        MPI_Gather(local_a, n/size, MPI_DOUBLE, array, n/size, MPI_DOUBLE, root, comm);
        time = MPI_Wtime() - time;
        commT += time;

        printf("Processor %d commT = %lf compT = %lf\n", rank, commT, compT);
        return MPI_SUCCESS;

    }
int MPI_Exchange( int n, double * array, int rank1, int rank2, MPI_Comm comm ) {
  int rank, size, result, i, tag1 = 0, tag2 = 1;
  double * b = ( double * ) calloc( n, sizeof( double ) );
  double * c;
  MPI_Status status;
  MPI_Comm_rank( comm, &rank );
  MPI_Comm_size( comm, &size );
  //L8.6
  if( rank == rank1 )
  {
    result = MPI_Send( &array[ 0 ], n, MPI_DOUBLE, rank2, tag1, comm );
    result = MPI_Recv( &b[ 0 ], n, MPI_DOUBLE, rank2, tag2, comm, &status );
    c = merge( n, array, n, b );
    for( i = 0; i < n; i++ )
    {
      array[ i ] = c[ i ];
    }
  }
  else if( rank == rank2 )
  {
    result = MPI_Recv( &b[ 0 ], n, MPI_DOUBLE, rank1, tag1, comm, &status );
    result = MPI_Send( &array[ 0 ], n, MPI_DOUBLE, rank1, tag2, comm) ;
    c = merge( n, array, n, b );
    for( i =0; i < n; i++ )
    {
      array[ i ] = c[ i + n ];
    }
  }
  return MPI_SUCCESS;
}
int MPI_Is_sorted(int n, double * array, int * answer, int root, MPI_Comm comm){

        // find rank and size
        int rank, size;
        MPI_Comm_rank( comm, &rank );
        MPI_Comm_size( comm, &size );

        // allocate space for size elements in first and last
        double * first = (double *)calloc(size, sizeof(double));
        double * last  = (double *)calloc(size, sizeof(double));

        // gather first and last
        MPI_Gather(&array[0],  1, MPI_DOUBLE, first, 1, MPI_DOUBLE, root, comm);
        MPI_Gather(&array[-1], 1, MPI_DOUBLE, last , 1, MPI_DOUBLE, root, comm);

        // if proc root then test
        if(rank == root){
            *answer = 1;
            for(int i=1;i<size;i++){
                if(last[i-1] > first[i]){
                    *answer = 0; break;
                }
            }
        }

        //
        MPI_Bcast(answer, 1, MPI_INT, root, comm);

        return MPI_SUCCESS;

    }


//notes
double * merge( int n, double * a, int m, double * b ) {
   int i, j, k;
   double * c = ( double * ) calloc( n + m, sizeof( double ) );
   for( i=j=k=0; ( i < n ) && ( j < m ); )
   {
      if( a[ i ] <= b[ j ] )
      {
        c[ k++ ] = a[ i++ ];
      }
      else
      {
        c[ k++ ] = b[ j++ ];
      }
   }
  if( i == n )
  {
    for( ; j < m; )
    {
      c[ k++ ] = b[ j++ ];
    }
  }
  else
  {
    for( ; i < n; )
    {
      c[ k++ ] = a[ i++ ];
    }
  }
  return c;
}
//notes
void merge_sort( int n, double * a ) {
  double * c;
  int i;
  if ( n <= 1 )
  {
    return;
  }
  if( n == 2 )
  {
    if( a[ 0 ] > a[ 1 ] )
    {
      swap( &a[ 0 ], &a[ 1 ] );
    }
    return;
  }
  merge_sort( n / 2, a );
  merge_sort( n - n / 2, a + n / 2 );
  c = merge( n / 2, a, n - n / 2, a + n / 2);
  for( i = 0; i < n; i++ )
  {
    a[ i ] = c[ i ];
  }
}
//notes
void swap ( double * a, double * b ) {
   double temp;
   temp = *a;
   *a = *b;
   *b = temp;
}