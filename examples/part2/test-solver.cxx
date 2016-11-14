#include <iostream>
#include <limits>
#include <cassert>
#include <numcxx/numcxx.hxx>


int main(void)
{
    auto pA=numcxx::TMatrix<double>::create(
        {   {2,3,5},
            {3,2,3},
            {9,5,7}});
    
    auto pLU=numcxx::TSolverLapackLU<double>::create(pA);
    pLU->update();
    auto pG=numcxx::TArray1<double>::create({1,2,3});
    auto pU=pG->clone();
    auto pF=pG->clone();

    auto &F=*pF;
    auto &G=*pG;
    auto &U=*pU;

    pLU->solve(U,G);
    pA->apply(U,F);
    double residual=normi(F-G);

    std::cout << "residual:" << residual << std::endl;
    assert(residual<1.0e4*std::numeric_limits<double>::epsilon());

    std::cout << U <<  std::endl;
}
