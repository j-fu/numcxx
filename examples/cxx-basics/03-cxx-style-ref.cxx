/// \example 03-cxx-style-ref.cxx
///
/// Working with C++ STL vectors passed as references
///
/// The class std::vector contains size information and can be
/// instantiated for any data type. We pass a reference to the vector
/// to the functions using it.
///
#include <cstdio>
#include <vector>


void initialize(std::vector<double>& x)
{
  const int n=x.size();
  for (int i=0;i<n;i++) x[i]= 1.0/(double)(1+n-i);
}

double sum_elements(std::vector<double>& x)
{
  double sum=0;
  for (int i=0;i<x.size();i++)sum+=x[i];
  return sum;
}

int main()
{
  const int n=12345678;
  
  // Instantiate a vector of doubles of given size.
  // The class object itself is placed on the stack,
  // but the data is placed on the heap.
  std::vector<double> x(n);

  initialize(x);
  double s=sum_elements(x);
  printf("sum=%e\n",s);
  
  // When leaving this scope, x is destroyed, and with
  // this also the memory on the heap is freed.
  return 0;
}
