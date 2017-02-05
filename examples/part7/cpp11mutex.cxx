#include <iostream>
#include <thread>
#include <mutex>

std::mutex mtx;

void call_from_thread(int tid) {
  mtx.lock();
  std::cout << "Launched by thread " << tid << std::endl;
  mtx.unlock();
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
  
  mtx.lock();
  std::cout << "Launched from the main\n";
  mtx.unlock();
  
  //Join the threads with the main thread
  for (int i = 0; i < num_threads; ++i) 
  {
    t[i].join();
  }
  
  return 0;
}

