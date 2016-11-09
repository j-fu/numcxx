#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include "numcxx/numcxx.hxx"
#include "numcxx/util.hxx"


void xsolve(numcxx::TLinSolver<double> &solver, numcxx::TArray<double>& u,const numcxx::TArray<double>& f)
{
   solver.solve(u,f);
}

int main()
{

    auto pA=numcxx::TMatrix<double>::create(
        {   {2,3,5},
            {3,2,3},
            {9,5,7}});
    
    auto pLU=std::make_shared<numcxx::TSolverLapackLU<double>>(pA);
    pLU->update();
    auto pG=numcxx::TArray1<double>::create({1,2,3});
    auto pU=pG->clone();
    xsolve(*pLU,*pU,*pG);
    std::cout<< *pU << std::endl;

    auto pInv=pA->calculate_inverse();
    std::cout << *pInv << std::endl;
    // auto pX=numcxx::linspace<double>(0,1,101);
    // std::cout << *pX;

}

