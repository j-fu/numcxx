/// \example 01-c-style-stack.cxx
///
/// Working with "plain" c-style arrays placed on the stack
///
/// Demonstrate  the effect of placing a large array on the stack.
/// 
/// Functions using plain c style arrays need to have
/// both the pointer to the data and the size information
/// as parameters.
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
  // Here, we place the array on the stack
  const int n=12345678;
  double x[n];

  initialize(x,n);
  double s=sum_elements(x,n);
  printf("sum=%e\n",s);

  // Before returning,  the space on the stack is freed
  return 0;
}
