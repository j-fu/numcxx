#include <iostream>
#include <thread>

//This function will be called from a thread

void call_from_thread(int tid) {
    std::cout << "Launched by thread " << tid << std::endl;
}

int main (int argc, char *argv[])
{
  
  int num_threads=1;
  if (argc>1) num_threads=atoi(argv[1]);
  
  std::thread t[num_threads];
  
  //Launch a group of threads
  for (int i = 0; i < num_threads; ++i) 
  {
    t[i] = std::thread(call_from_thread, i);
  }
  
  std::cout << "Launched from the main\n";
  
  //Join the threads with the main thread
  for (int i = 0; i < num_threads; ++i) 
  {
    t[i].join();
  }
  
  return 0;
}

