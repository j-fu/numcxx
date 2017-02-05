#include <iostream>
#include <cstdlib>
#include <omp.h>

void call_from_thread() {
  int tid=omp_get_thread_num();
#pragma omp critical
  {
    std::cout << "Launched by thread " << tid << std::endl;
  }
}

int main (int argc, char *argv[])
{
  
  int num_threads=1;
  if (argc>1) num_threads=atoi(argv[1]);
  
  
#pragma omp parallel for
  for (int i = 0; i < num_threads; ++i) 
  {
    call_from_thread();
  }
  return 0;
}


