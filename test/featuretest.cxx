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
    assert(min(*a)==3);
    assert(max(*a)==21);
    assert(sum(*a)==120);
    assert(norm1(*a)==120);
    assert(normi(*a)==21);
    (*a)-=4;
    (*a)(0)+=4;
    (*a)(2)+=1;
    (*a)(3)+=3;
    assert(norm2(*a)==32); //exact value, power of 2, works as well for integer types

    auto b=a->clone();
    *b=*a;
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
        }
}

template <typename T > void matrixtest(void)
{

    auto A=numcxx::TMatrix<double>::create(
        {   {2,3,5},
            {3,2,3},
            {9,5,7}});
    
    auto LU=std::make_shared<numcxx::TSolverLapackLU<double>>(A);
    LU->update();
    auto g=numcxx::TArray1<double>::create({1,2,3});
    auto u=g->clone();
    auto f=g->clone();
    LU->solve(*u,*g);
    A->apply(*u,*f);
    assert(normi(*f-(*g))<1.0e4*std::numeric_limits<T>::epsilon());
}


template <typename T > void expressiontest(void)
{
    const int N=3;
    numcxx::TArray1<T> A(N);
    numcxx::TArray1<T> B(N);
    numcxx::TArray1<T> C(N);
    numcxx::TArray1<T> D(N);

    numcxx::TMatrix<T> Mat{{1,2,3},{4,5,6},{7,8,9}};
    auto X=Mat.clone();
    auto &x=*X;
    A=0.0;

    x=Mat;

    B=1;
    C=21;
    D=3;
    A=B+5;
    A=6+B;
    D=A*4;
    A=3*(B+C/3)-4*D/2;
    assert(norm1(A)==72);



    B=3*(Mat*A);
    assert(norm1(B)==3240);
    Mat.apply(A,C);
    B=Mat*A;
    assert(norm1(B)==norm1(C));
    

    // std::vector<T> E{33,44,55};
    // B=A+E;
    // B+=Mat*E;
    // assert(norm1(B)==2106);
    // B=E;
    // assert(norm1(B)==132);
}

int main()
{

    arraytest<double>();
    arraytest<long double>();
    arraytest<float>();
    arraytest<int>();
    arraytest<long int>();
    arraytest<short int>();

    matrixtest<double>();
    matrixtest<float>();

    expressiontest<double>();
    expressiontest<long double>();
    expressiontest<float>();
    expressiontest<int>();
    expressiontest<long int>();
    expressiontest<short int>();


}
