#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include "numcxx/numcxx.hxx"


void xsolve(numcxx::TLinSolver<double> &solver, numcxx::TArray<double>& u,const numcxx::TArray<double>& f)
{
   solver.solve(u,f);
}

int main()
{

    auto A=numcxx::TMatrix<double>::create(
        {   {2,3,5},
            {3,2,3},
            {9,5,7}});
    
    auto LU=std::make_shared<numcxx::TSolverLapackLU<double>>(A);
    LU->update();
    auto g=numcxx::TArray1<double>::create({1,2,3});
    auto u=g->clone();
    xsolve(*LU,*u,*g);
    std::cout<< u << std::endl;
}

