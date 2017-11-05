/// \example 02-c-style-heap.cxx
///
/// Working with "plain" c-style arrays placed on the heap
///

#include <cstdio>


void initialize(double *x, int n)
{
  for (int i=0;i<n;i++) x[i]= 1.0/(double)(1+n-i);
}

double sum_elements(double *x, int n)
{
  double sum=0;
  for (int i=0;i<n;i++) sum+=x[i];
  return sum;
}

int main()
{
  const int n=12345678;

  // Allocate an array on the heap
  double* x=new double[n];

  initialize(x,n);
  double s=sum_elements(x,n);
  printf("sum=%e\n",s);
  
  // Free the space previously allocated.
  delete[] x;
  return 0;
}

