#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include "numcxx/numcxx.hxx"
#include "netlib/netlib.hxx"

void tdense(void)
{

    auto A=numcxx::TMatrix<double>::create(
        {   {2,3,5},
            {3,2,3},
            {9,5,7}});

    auto pSolver=std::make_shared<numcxx::TSolverLapackLU<double>>(A);
    pSolver->update();
    auto pG=numcxx::TArray1<double>::create({1,2,3});
    auto pU=pG->clone();
    pSolver->solve(*pU,*pG);
    std::cout<< *pU << std::endl;

}

void tsparse()
{

    auto pA=numcxx::TSparseMatrix<double>::create(
        {   {2,3,5},
            {3,2,3},
            {9,5,7}});

    auto pDA=pA->copy_as_dense();

    std::cout << *pDA << std::endl;

    auto pSolver=std::make_shared<numcxx::TSolverUMFPACK<double>>(pA);
    pSolver->update();
    auto pG=numcxx::TArray1<double>::create({1,2,3});
    auto pU=pG->clone();
    pSolver->solve(*pU,*pG);
    std::cout<< *pU << std::endl;



    int n=10;

    auto pM=std::make_shared<numcxx::DSparseMatrix>(n);
    auto &M=*pM;
    auto F=numcxx::DArray1(n);
    F=1.0;
    auto U=F;
    auto V=F;



    for (int i=0;i<n;i++)
    {
        M(i,i)=3.0;
        if (i>0)
            M(i,i-1)=-1;
        if (i<n-1)
            M(i,i+1)=-1;
    }


    auto pMInv=M.calculate_inverse();
    std::cout << *pMInv << std::endl;
    

    auto pDM=M.copy_as_dense();

    double tic,toc;
    tic=numcxx::cpu_clock();
    auto pLapack=numcxx::DSolverLapackLU::create(pDM);
    pLapack->solve(U,F);
    toc=numcxx::cpu_clock();
    std::cout << "seconds for dense:" << toc-tic << std::endl;
    
    tic=numcxx::cpu_clock();
    auto pUmfpack=std::make_shared<numcxx::DSolverUMFPACK>(pM);
    pUmfpack->solve(V,F);
    toc=numcxx::cpu_clock();
    std::cout << "seconds for sparse:" << toc-tic << std::endl;
    std::cout << "difference sparse-dense: "<< normi(U-V) << std::endl;
    std::cout << U << std::endl;
}



void titer()
{
    int n=100;

    auto pM=std::make_shared<numcxx::DSparseMatrix>(n);
    auto &M=*pM;
    auto F=numcxx::DArray1(n);
    F=1.0;
    auto U=F;

    for (int i=0;i<n;i++)
    {
        M(i,i)=3.0;
        if (i>0)
            M(i,i-1)=-1;
        if (i<n-1)
            M(i,i+1)=-1;
    }
    M.flush();
    auto pJacobi=std::make_shared<numcxx::TPreconJacobi<double>>(pM);

    auto R=U;
    auto V=U;
    int max_iter=1000;
    for (int niter=0;niter<max_iter;niter++)
    {
        R=M*U-F;
        pJacobi->solve(V,R);
        U-=V;
    }
    std::cout << U << std::endl;

    V=F;

    double tol=1.0e-10;
    netlib::BiCGSTAB(M,V,F,*pJacobi,max_iter,tol);
    std::cout << V << std::endl;
    
}

//CG(const Matrix &A, Vector &x, const Vector &b,
//   const Preconditioner &M, int &max_iter, Real &tol)


int main(void)
{
    int n=4;
    numcxx::TSparseMatrix<double> A(n);
    numcxx::TArray1<double>X(n);
    numcxx::TArray1<double>Y(n);
    numcxx::TArray1<double>Z(n);


    X=2.0;

    A(0,0)=10.0;
    A(1,1)=1.0;
    A(1,2)=100.0;
    A(2,2)=2;
    A(3,3)=3;
    A.flush();
    A.apply(X,Y);
    std::cout <<Y << std::endl;

    auto pM=A.copy_as_dense();
    pM->apply(X,Z);
    std::cout <<*pM << std::endl;
    std::cout <<Z << std::endl;

    tdense();
    tsparse();
//    titer();
//    std::cout << x->Matrix.shape(0) << std::endl;
}
