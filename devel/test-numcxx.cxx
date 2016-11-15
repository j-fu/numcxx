#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include "numcxx/numcxx.hxx"
int test_numcxx2(void);

void speedtest(int n)
{
  int nrun=200;
  auto XA=numcxx::TArray1<double>::create(n);
  auto &xa=*XA;
  double *ca=new double[n];
  ca[0]=1;
  xa(0)=1;
        
  double tc = (double)std::clock()/(double)CLOCKS_PER_SEC;
  for (int irun=0;irun<nrun;irun++)
    for (int i=1;i<n;i++)
      ca[i]=(double)(i*irun)+ca[i-1];
  double dc = (double)std::clock()/(double)CLOCKS_PER_SEC-tc;
  printf("dc=%e s\n",dc);

  
  double tx = (double)std::clock()/(double)CLOCKS_PER_SEC;
  for (int irun=0;irun<nrun;irun++)
    for (int i=1;i<n;i++)
      {
        xa(i)=(double)(i*irun)+xa(i-1);
      }
  double dx = (double)std::clock()/(double)CLOCKS_PER_SEC-tx;


  printf("dx=%e s\n",dx);
  
  delete[] ca;
}

template <typename T> 
void xtest(numcxx::TArray<T> & x)
{
     std::cout << x.size() <<  std::endl;
}

int main (int argc, const char *argv[])
{
    const int N=4;

    auto v=std::make_shared<std::vector<double>>(10);
    for (int i=1;i<10;i++) (*v)[i]=-i;

    auto xv=std::make_shared<numcxx::DArray1>(v);
    std::cout << *xv << std::endl;
    xtest(*xv);

    auto A=numcxx::DArray1::create({3,4,5});

    auto B=A->clone();
    auto C=A->clone();
    numcxx::DArray1::operate([](double &a, double &b){b=std::sin(a);},*A,*B);
    
    numcxx::DArray1::operate([](double &a, double &b, double &c){c=2.0*a+b;},*A,*B,*C);

    std::cout << "A:\n"<<A << std::endl;

    std::cout << "B:\n"<<B << std::endl;

    std::cout << "C:\n"<<C<< std::endl;

    auto M=numcxx::DMatrix::create(
        {
            {3,4,7},
            {5,6,-1},
            {0,9,18}
        });
    std::cout << "M:\n"<<*M<< std::endl;
    auto a=numcxx::DArray1::create({1,2,3});
    auto b=a->clone();
    auto solver=numcxx::TSolverLapackLU<double>::create(M);
    solver->update();
    solver->solve(*b,*a);

    std::cout<< "b:\n"<< b << std::endl;

    auto c=b->clone();
    M->apply(*b,*c);
    
    std::cout<< "c:\n"<< c << std::endl;
    
 
    std::cout<< test_numcxx2() << std::endl;
    speedtest(2000000);
}


