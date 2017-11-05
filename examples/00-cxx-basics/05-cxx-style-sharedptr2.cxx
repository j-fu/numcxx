///
/// \example 05-cxx-style-sharedptr2.cxx
///
/// Working with C++ vectors allocated as smart pointers
/// but passed as references.  Demonstrate as well the use of 
/// [std::make_shared](http://www.cplusplus.com/reference/memory/make_shared/).
/// 
#include <cstdio>
#include <vector>
#include <memory>


/// We pass (as in example 03) a reference to
/// the vector to the functions using it. This is possible as long as the
/// function itself does not store the vector in another class.

void initialize(std::vector<double> &x)
{
  const int n=x.size();
  for (int i=0;i<n;i++) x[i]= 1.0/(double)(1+n-i);
}

double sum_elements(std::vector<double> & x)
{
  double sum=0;
  for (int i=0;i<x.size();i++)sum+=x[i];
  return sum;
}

int main()
{
  const int n=12345678;

  // This is a shortcut to the cumbersome statement from
  // Example 04. We have to write the data type only once,
  // as a template argument to std::make_shared<>.
  auto x=std::make_shared<std::vector<double>>(n);
  
  // Dereference the smart pointer here in order to 
  // obtain a reference which can be used as an argumet to the functions.
  initialize(*x);
  double s=sum_elements(*x);
  printf("sum=%e\n",s);

  // The shared pointer is destroyed, reference counter decremented
  // and, as it becomes 0, the allocated memory is released.
  return 0;
}

