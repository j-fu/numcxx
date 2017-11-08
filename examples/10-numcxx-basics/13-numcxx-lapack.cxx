/// 
///  \example 13-numcxx-lapack.cxx
///
/// Use of lapack solver for dense matrices
/// 
#include <iostream>
#include <limits>
#include <numcxx/numcxx.hxx>


int main(void)
{
  int n=5000;
  numcxx::DMatrix A(n,n);
  numcxx::DArray1 F(n);
  numcxx::DArray1 U(n);
  
  
  F=1.0;
  for (int i=0;i<n;i++)
  {
    A(i,i)=3.0;
    if (i>0) A(i,i-1)=-1;
    if (i<n-1) A(i,i+1)=-1;
  }
  
  numcxx::DSolverLapackLU LapackSolver(A);
  LapackSolver.solve(U,F);
  double residual=normi(A*U-F);
  
  std::cout << "residual:" << residual << std::endl;
}

