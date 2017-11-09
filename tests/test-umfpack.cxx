#include <iostream>
#include <limits>
#include <cassert>
#include <numcxx/numcxx.hxx>


int main(void)
{
    int n=5000;
    auto pM=numcxx::DSparseMatrix::create(n,n);
    auto pF=numcxx::DArray1::create(n);
    auto pU=numcxx::DArray1::create(n);
    
    auto &M=*pM;
    auto &F=*pF;
    auto &U=*pU;
    
    F=1.0;
    for (int i=0;i<n;i++)
    {
        M(i,i)=3.0;
        if (i>0) M(i,i-1)=-1;
        if (i<n-1) M(i,i+1)=-1;
    }
    pM->flush();
    
    auto pUmfpack=numcxx::DSolverUMFPACK::create(pM);
    pUmfpack->update();
    pUmfpack->solve(U,F);
    double residual=normi(M*U-F);
    
    std::cout << "residual:" << residual << std::endl;
    assert(residual<1.0e-14);
}


