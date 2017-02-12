# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <mpi.h>
# include <vector>
#include <chrono>
#include <thread>


using namespace std;

int main ( int argc, char *argv[] )
{
  const unsigned long N=2000001001;
  int iproc;
  int nproc;
  double tstart;
  MPI_Status status;
  unsigned long sum_proc;
  unsigned long sum_all;
  const int master=0;
  const int tag=1;

  //  Initialize MPI.
  MPI_Init ( &argc, &argv );
  
  //  Get the number of processes.
  MPI_Comm_size ( MPI_COMM_WORLD, &nproc );
  
  // Create index vector for processes
  std::vector<unsigned long> idx(nproc+1);
  
  //  Determine the rank (number) of this process.
  MPI_Comm_rank ( MPI_COMM_WORLD, &iproc );
  
  

  if ( iproc == 0 )
  {
    tstart=MPI_Wtime();
    cout << "The number of processes available is " << nproc << "\n";
  }
  
  // The master process initializes the numbers to be handeld
  unsigned long chunksize=N/nproc;
  if ( iproc == 0 ) 
  {
    unsigned long chunkstart=0;
    for (int i=0;i<nproc;i++)
    {
      idx[i]=chunkstart;
      chunkstart+=chunksize;
    } 
    idx[nproc]=N;
  }

  //  The master process broadcasts the computed indices
  //  to all the other processes.
  MPI_Bcast ( idx.data(), nproc+1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
  // From here on, all processes have the same data in the vector idx.
  
  // Initialize and calculate the partial sum corresponding
  // to process number
  sum_proc=0;
  for (unsigned long i=idx[iproc];i<idx[iproc+1];i++)
    sum_proc+=i;

  cout << "Proc  " << iproc<< ": sum( " << idx[iproc] <<  " ... " << idx[iproc+1]-1 <<" ) = " << sum_proc << "\n";

  //  Each worker process sends its sum back to the master process.
  if ( iproc != 0  ) 
  {
    MPI_Send ( &sum_proc, 1, MPI_UNSIGNED_LONG, master, tag, MPI_COMM_WORLD );
  }
  else 
  {
    sum_all = sum_proc;
    // Master process receives data from the others, sequence is undefined
    for ( int i = 1; i < nproc; i++ ) 
    {
      MPI_Recv ( &sum_proc, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, tag,  MPI_COMM_WORLD, &status );
      cout << "Received result "<< sum_proc << " from proc "<< status.MPI_SOURCE << endl;
      sum_all = sum_all + sum_proc;
    }
  }
  double tend= MPI_Wtime();
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  if ( iproc == 0 )
  {
    cout << "         sum( " << 0 <<  " ... " << N-1 <<" ) = " << sum_all << endl;
    cout << "        check:            "<< (N*(N-1))/2<<"\n";
    cout << "wall clock time/s: " << tend-tstart; 
  }

  
  //  Terminate MPI.
  MPI_Finalize ( );

  if ( iproc == 0 )
  {
    cout << "  Normal end of execution.\n";
  }
  return 0;


}
