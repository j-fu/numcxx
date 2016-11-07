#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include "numcxx/numcxx.hxx"
#include "numcxx/tsparsematrix.hxx"
#include "numcxx/tsolver-umfpack.hxx"

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

void tsparse(void)
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

}



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
}
