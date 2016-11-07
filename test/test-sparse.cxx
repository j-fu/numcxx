#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include "numcxx/numcxx.hxx"
#include "numcxx/tsparsematrix.hxx"


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

  
}
