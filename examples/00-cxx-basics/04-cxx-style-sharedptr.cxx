/// \example 04-cxx-style-sharedptr.cxx
///
/// Working with C++ vectors passed via shared pointers
///
/// This is one way to work with smart pointers of the type
/// [std::shared_ptr](http://www.cplusplus.com/reference/memory/shared_ptr/)
///
#include <cstdio>
#include <vector>
#include <memory>

// The smart pointer is passed by value, it contains
// a pointer to the object.
void initialize(std::shared_ptr<std::vector<double>> x)
{
  // Here we need to dereference the smart pointer
  // in order to access the data. This looks a  bit
  // confusing...
  const int n=x->size();
  for (int i=0;i<n;i++) (*x)[i]= 1.0/(double)(1+n-i);
}

double sum_elements(std::shared_ptr<std::vector<double>> x)
{
  double sum=0;
  for (int i=0;i<x->size();i++)sum+=(*x)[i];
  return sum;
}

int main()
{
  const int n=12345678;

  // This declares  and creates a shared pointer from an object
  // allocated with new. Quite cumbersome...
  std::shared_ptr<std::vector<double>> x=std::shared_ptr<std::vector<double>>(new std::vector<double>(n));
  initialize(x);
  double s=sum_elements(x);
  printf("sum=%e\n",s);

  // The shared pointer is destroyed, reference counter decremented
  // and, as it becomes 0, the allocated memory is released.
  return 0;
}
