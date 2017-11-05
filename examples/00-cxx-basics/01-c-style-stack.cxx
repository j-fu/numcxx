/**
   \example 01-c-style-stack.cxx 

   Working with "plain" c-style arrays placed on the stack.
*/

#include <cstdio>

// Initialize the vector with some data
//
// Functions using plain c style arrays need to have
// both the pointer to the data and the size information
// as parameters.
void initialize(double *x, int n)
{
  for (int i=0;i<n;i++) x[i]= 1.0/(double)(1+n-i);
}

// Sum up the elements of the vector and return the value
double sum_elements(double *x, int n)
{
  double sum=0;
  for (int i=0;i<n;i++) sum+=x[i];
  return sum;
}

int main()
{
  const int n=12345678;
  //  Declare the array and place it on the stack
  double x[n];

  initialize(x,n);
  double s=sum_elements(x,n);
  printf("sum=%e\n",s);

  // Before returning,  the space on the stack is freed
  return 0;
}
