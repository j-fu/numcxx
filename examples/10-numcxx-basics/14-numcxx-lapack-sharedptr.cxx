/// 
///  \example 14-numcxx-lapack-sharedptr.cxx
///
/// Use of lapack solver for sparse matrices, shared pointer idiom
/// 

#include <iostream>
#include <limits>
#include <numcxx/numcxx.hxx>


int main(void)
{
    int n=5000;
    auto pA=numcxx::DMatrix::create(n,n);
    auto pF=numcxx::DArray1::create(n);
    auto pU=numcxx::DArray1::create(n);
    
    // Here, we need to derefernce the smart pointer to references
    // in order to use the () operator in a convenient way.
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
    
    auto pLapack=numcxx::DSolverLapackLU::create(pA);
    pLapack->solve(U,F);
    double residual=normi(A*U-F);
    
    std::cout << "residual:" << residual << std::endl;
}

