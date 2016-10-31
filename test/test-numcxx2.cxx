#include "numcxx/numcxx.hxx"


int test_numcxx2(void)
{
  auto a=numcxx::TArray1<double>::create(10);

  return a->size();
}

