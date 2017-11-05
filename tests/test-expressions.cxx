#include <iostream>
#include <numcxx/numcxx.hxx>
#include <stdexcept> 

//
// Check all basic operations for expression templates
//

template <typename T1, typename T2, typename T3>
void check(T1& C1, T2& C2, T3& C3)
{
  int n=C1.size();
  std::cout << typeid(T1).name() << std::endl;
  std::cout << typeid(T2).name() << std::endl;
  std::cout << typeid(T3).name() << std::endl;
  
  for (int i=0;i<n;i++)
  {
    //std::cout << i << " " << C1[i] <<" " << C2[i] <<" " << C3[i] << std::endl;
    if (C1[i]!=C2[i]) 
    {
      std::cout <<"error in expression template" << std::endl;
      exit(1);
    }
    if (C1[i]!=C3[i]) 
    {
      std::cout <<"error in expression template" << std::endl;
      exit(1);
    }
  }
  std::cout << std::endl;
}

int main()
{
    const int n=100;
    
    numcxx::DArray1 A(n);
    numcxx::DArray1 B(n);
    numcxx::DArray1 C1(n);
    for (int i=0;i<n;i++)
    {
      A(i)=10+i;
      B(i)=10+n-i-1;
    }
    
    
    {
      std::cout << "operation A+B" << std::endl;
      for (int i=0;i<n;i++) 
        C1(i)=A(i)+B(i);
      numcxx::DArray1 C2=A+B;
      auto C3=A+B;
      check(C1,C2,C3);
    }

    {
      std::cout << "operation A-B" << std::endl;
      for (int i=0;i<n;i++) 
        C1(i)=A(i)-B(i);
      numcxx::DArray1 C2=A-B;
      auto C3=A-B;
      check(C1,C2,C3);
    }

    {
      double a=3.0;
      std::cout << "operation a*B" << std::endl;
      for (int i=0;i<n;i++) 
        C1(i)=a*B(i);
      numcxx::DArray1 C2=a*B;
      auto C3=3.0*B;
      check(C1,C2,C3);
      C1=C3*2.0; // this triggers optimizing away scalar reference in g++
    }

    {
      double b=3.0;
      std::cout << "operation A*b" << std::endl;
      for (int i=0;i<n;i++) 
        C1(i)=A(i)*b;
      numcxx::DArray1 C2=A*b;
      auto C3=A*b;
      check(C1,C2,C3);
      C1=C3*2.0; // this triggers optimizing away scalar reference in g++
    }


    {
      double a=3.0;
      std::cout << "operation a+B" << std::endl;
      for (int i=0;i<n;i++) 
        C1(i)=a+B(i);
      numcxx::DArray1 C2=a+B;
      auto C3=3.0+B;
      check(C1,C2,C3);
      C1=C3*2.0; // this triggers optimizing away scalar reference in g++
    }

    {
      double b=3.0;
      std::cout << "operation A+b" << std::endl;
      for (int i=0;i<n;i++) 
        C1(i)=A(i)+b;
      numcxx::DArray1 C2=A+b;
      auto C3=A+b;
      check(C1,C2,C3);
      C1=C3*2.0; // this triggers optimizing away scalar reference in g++
    }

    {
      double a=3.0;
      std::cout << "operation a-B" << std::endl;
      for (int i=0;i<n;i++) 
        C1(i)=a-B(i);
      numcxx::DArray1 C2=a-B;
      auto C3=3.0-B;
      check(C1,C2,C3);
      C1=C3*2.0; // this triggers optimizing away scalar reference in g++
    }

    {
      double b=3.0;
      std::cout << "operation A-b" << std::endl;
      for (int i=0;i<n;i++) 
        C1(i)=A(i)-b;
      numcxx::DArray1 C2=A-b;
      auto C3=A-b;
      check(C1,C2,C3);
      C1=C3*2.0; // this triggers optimizing away scalar reference in g++
    }

    {
      double b=3.0;
      std::cout << "operation A/b" << std::endl;
      for (int i=0;i<n;i++) 
        C1(i)=A(i)/b;
      numcxx::DArray1 C2=A/b;
      auto C3=A/b;
      check(C1,C2,C3);
      C1=C3*2.0; // this triggers optimizing away scalar reference in g++
    }


    {
      std::cout << "operation M*A" << std::endl;
      numcxx::DMatrix M(n,n);
      for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
          M(i,j)=10*i+3*j+1;
      M.apply(A,C1);
      numcxx::DArray1 C2=M*A;
      auto C3=M*A;
      check(C1,C2,C3);
    }

    {
      
      std::cout << "operation M*A (sparse)" << std::endl;
      numcxx::DSparseMatrix M(n,n);
      for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
          M(i,j)=10*i+3*j+1;
      M.flush();
      M.apply(A,C1);
      numcxx::DArray1 C2=M*A;
      auto C3=M*A;
      check(C1,C2,C3);
    }


}

