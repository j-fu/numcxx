//
// /net/wir/numcxx/examples/20-graphics/20-numcxx-matplotlib.py
//
#include <numcxx/numcxx.hxx>
#include <cstdio>
#include <cmath>


int main(void)
{
  const double x0=0.0;
  const double x1=1.0;
  double h=0.01;
  double t0=0.0;
  double t1=10.0;
  double dt=0.01;

  const int N=1+ceil((x1-x0)/h);
  h=(x1-x0)/(double)(N-1);

  numcxx::DArray1 X(N);
  numcxx::DArray1 U(N);

  X(0)=x0;
  for (int i=1;i<N;i++)
  {
    X(i)=X(i-1)+h;
  }
  
  for (int i=0; i<N; i++)
  {
    U(i)=sin(20.0*X(i));
  }
  numcxx::savetxt("20-plot-X.dat",X);
  numcxx::savetxt("20-plot-U.dat",U);
  
}


