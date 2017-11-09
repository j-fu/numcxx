/// 
///  \example 15-numcxx-umfpack.cxx
///
/// Use of umfpack solver for sparse matrices
/// 


#include <iostream>
#include <limits>
#include <cassert>
#include <numcxx/numcxx.hxx>


int main(void)
{
    int n=5000;
    numcxx::DSparseMatrix A(n,n);
    numcxx::DArray1 U(n);
    numcxx::DArray1 F(n);
    
    F=1.0;
    for (int i=0;i<n;i++)
    {
        A(i,i)=3.0;
        if (i>0) A(i,i-1)=-1;
        if (i<n-1) A(i,i+1)=-1;
    }
    A.flush();
    
    numcxx::DSolverUMFPACK UmfpackSolver(A);
    UmfpackSolver.update(A);
    UmfpackSolver.solve(U,F);
    double residual=normi(A*U-F);
    
    std::cout << "residual:" << residual << std::endl;
}


