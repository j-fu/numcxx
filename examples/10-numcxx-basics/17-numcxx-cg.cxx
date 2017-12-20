/// 
///  \example 17-numcxx-cg.ccxx
///
/// Use of cg for sparse matrices with shared ptr idiom
/// 

#include <iostream>
#include <limits>
#include <cassert>
#include <numcxx/numcxx.hxx>
#include "netlib/netlib.hxx"


int main(void)
{
  int n=5000; // problem size
  int max_iter=100; // maximum number of iterations
  double tol=1.0e-14; // tolerance
  
  auto pA=numcxx::DSparseMatrix::create(n,n);
  auto pF=numcxx::DArray1::create(n);
  auto pU=numcxx::DArray1::create(n);
  
  auto &A=*pA;
  auto &F=*pF;
  auto &U=*pU;
  
  F=1.0;
  for (int i=0;i<n;i++)
  {
    A(i,i)=3.0;
    if (i>0) A(i,i-1)=-1;
    if (i<n-1) A(i,i+1)=-1;
  }
  pA->flush();
  
  
  auto pJacobi=numcxx::DPreconJacobi::create(pA);
  numcxx::TPreconJacobi<double> JacobiSolver(pA);

  
  U=0.0;
  netlib::CG(A,U,F,*pJacobi,max_iter,tol);
  double residual=normi(A*U-F);
  std::cout << "residual:" << residual << std::endl;
}



