#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include "numcxx/numcxx.hxx"



int main()
{
    const int N=10;
    numcxx::DArray1 A(N);
    numcxx::DArray1 B(N);
    numcxx::DArray1 C(N);
    numcxx::DArray1 D(N);
        
    A=0.0;
    B=1.0;
    C=2.0;
    D=3.25;
    A=B;

    A=3*(B+C/3.0)-4.0*D/2.0;

    std::cout<< A;

}
