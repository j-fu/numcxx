#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cassert>
#include <limits>
#include "numcxx/numcxx.hxx"



template <typename T > void arraytest(void)
{
    auto a= numcxx::TArray1<T>::create({1,2,3,4,5,6,7,8,9,10});
    assert(a->ndim()==1);
    assert(a->size()==10);
    assert(a->shape(0)==10);
    for (int i=0;i<a->size();i++) assert((*a)(i)==i+1);
    for (int i=0;i<a->size();i++) assert((*a)[i]==i+1);
    
    (*a)+=1;
    for (int i=0;i<a->size();i++) assert((*a)[i]==i+2);
    (*a)*=10;
    for (int i=0;i<a->size();i++) assert((*a)[i]==10*(i+2));
    (*a)/=5;
    for (int i=0;i<a->size();i++) assert((*a)[i]==2*(i+2));
    (*a)-=1;
    for (int i=0;i<a->size();i++) assert((*a)[i]==2*(i+2)-1);
    assert(a->min()==3);
    assert(a->max()==21);
    assert(a->sum()==120);
    assert(a->norm1()==120);
    assert(a->normi()==21);
    (*a)-=4;
    (*a)(0)+=4;
    (*a)(2)+=1;
    (*a)(3)+=3;
    assert(a->norm2()==32); //exact value, power of 2, works as well for integer types

    auto b=a->clone();
    b->fill(*a);
    for (int i=0;i<a->size();i++) assert((*a)[i]==(*b)[i]);

    numcxx::TArray1<T>::operate([](T&a, T&b){a=2*b;},*a,*b);
    for (int i=0;i<a->size();i++) assert((*a)[i]==2*(*b)[i]);

    auto c=a->clone();
    numcxx::TArray1<T>::operate([](T&a, T&b, T&c){c=a+b;},*a,*b,*c);
    for (int i=0;i<a->size();i++) assert((*c)[i]=(*a)[i]+(*b)[i]);
    
    auto A= numcxx::TArray2<T>::create({{1,2,3},{4,5,6},{7,8,9}});
    for (int i=0;i<A->shape(0);i++)
        for (int j=0;j<A->shape(1);j++)
        {
            assert((*A)(i,j)==i*A->shape(1)+j+1);
            assert((*A)[i][j]==i*A->shape(1)+j+1);
        }
}

template <typename T > void matrixtest(void)
{
    auto M=numcxx::TMatrix<T>::create(
        {{1,1},
         {0,1}});
    auto f=numcxx::TArray1<T>::create({1,1});
    auto u=M->solve(*f);
    assert((*u)(0)==0);
    assert((*u)(1)==1);


    auto invM=M->inverse();
    assert((*invM)(0,0)==1);
    assert((*invM)(0,1)==-1);
    assert((*invM)(1,0)==0);
    assert((*invM)(1,1)==1);

    
    // this is an unimodular matrix with integer inverse, see
    // Weisstein, Eric W. "Unimodular Matrix." From MathWorld--A Wolfram Web Resource. 
    // http://mathworld.wolfram.com/UnimodularMatrix.html
    auto N=numcxx::TMatrix<T>::create(
        {   {2,3,5},
            {3,2,3},
            {9,5,7}});
    auto invN=N->inverse();

    auto g=numcxx::TArray1<T>::create({1,2,3});
    auto u1=N->solve(*g);
    auto u2=u1->clone();
    u2->fill(0);
    invN->apply(*g,*u2);
    (*u2)-=(*u1);
    
    assert(u2->normi()<50.0*std::numeric_limits<T>::epsilon());

    

}

int main()
{

    arraytest<double>();

    // arraytest<long double>();
    // arraytest<float>();
    // arraytest<int>();
    // arraytest<long int>();
    // arraytest<short int>();

    matrixtest<double>();
    matrixtest<float>();

}
