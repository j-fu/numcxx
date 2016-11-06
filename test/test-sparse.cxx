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
    numcxx::TArray1<double>x(n);
    auto y=x.clone();

    x=2.0;

    A(0,0)=10.0;
    A(1,1)=1.0;
    A(1,2)=100.0;
    A(2,2)=2;
    A(3,3)=3;
//    A.flush();

    A.apply(x,*y);
    std::cout <<*y << std::endl;
}
