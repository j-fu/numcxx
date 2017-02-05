/// Example with controlled  access
/// to variable via critical region - slow


#include <cstdio>
#include <numcxx/numcxx.hxx>
#include <omp.h>

int main(void)
{
  int N=1000000;
  auto pU=numcxx::DArray1::create(N);
  auto pV=numcxx::DArray1::create(N);

  auto &U=*pU;
  auto &V=*pV;
  U=10;
  V=20;

  for (int irun=1;irun<=1000; irun++)
  {
    double s=0.0;
#pragma omp parallel for
    for (int i=0;i<N;i++)
#pragma omp critical
    {
      s+=U(i)*V(i);
    }
    if (irun%100==0)
      printf("s=%g\n",s);
    
  }
}


