/// 
///  \example 16-numcxx-umfpack-sharedptr.cxx
///
/// Use of umfpack solver for sparse matrices with shared ptr idiom
/// 

#include <iostream>
#include <limits>
#include <cassert>
#include <numcxx/numcxx.hxx>


int main(void)
{
    int n=5000;
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
    
    auto pUmfpack=numcxx::DSolverUMFPACK::create(pA);
    pUmfpack->update();
    pUmfpack->solve(U,F);
    double residual=normi(A*U-F);
    
    std::cout << "residual:" << residual << std::endl;
}


