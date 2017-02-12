# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <mpi.h>
# include <vector>
#include <chrono>
#include <thread>
#include <unistd.h>



using namespace std;

int main ( int argc, char *argv[] )
{
  int iproc;
  int nproc;
  double tstart;
  MPI_Status status;
  
  //  Initialize MPI.
  MPI_Init ( &argc, &argv );
  
  //  Get the number of processes.
  MPI_Comm_size ( MPI_COMM_WORLD, &nproc );
  

  //  Determine the rank (number) of this process.
  MPI_Comm_rank ( MPI_COMM_WORLD, &iproc );

  if ( iproc == 0 )
  {
    tstart=MPI_Wtime();
    cout << "The number of processes available is " << nproc << "\n";
  }
 
  char hostname[80];
  gethostname(hostname,80);
  cout << "Hello from proc  " << iproc <<  " at host " << hostname << endl;

  double tend= MPI_Wtime();
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  if ( iproc == 0 )
  {
    cout << "wall clock time/s: " << tend-tstart << endl; 
  }

  //  Terminate MPI.
  MPI_Finalize ( );

  if ( iproc == 0 )
  {
    cout << "  Normal end of execution.\n";
  }
  return 0;
}
