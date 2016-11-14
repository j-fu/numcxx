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
    auto pR=numcxx::DArray1::create(n);
    auto pV=numcxx::DArray1::create(n);
    
    auto &M=*pM;
    auto &F=*pF;
    auto &U=*pU;
    auto &V=*pV;
    auto &R=*pR;


    F=1.0;
    for (int i=0;i<n;i++)
    {
        M(i,i)=3;
        if (i>0) M(i,i-1)=-1;
        if (i<n-1) M(i,i+1)=-1;
    }
    pM->flush();
    auto pJacobi=numcxx::DPreconJacobi::create(pM);
    pJacobi->update();
    double residual_norm=0.0;
    U=0.0;
    int niter=1000;
    for (int i=0;i<niter;i++)
    {
        R=M*U-F;
        residual_norm=normi(R);    
        if (residual_norm<1.0e-15) break;
        pJacobi->solve(V,R);
        U-=V;
    }
    
    std::cout << "residual:" << residual_norm << std::endl;

}

