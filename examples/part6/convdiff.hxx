#include <iostream>
#include <limits>
#include <cassert>
#include <numcxx/numcxx.hxx>
#include <cmath>

namespace convdiff
{

  inline double B(double x)
  {
    if (std::fabs(x)<1.0e-10) return 1.0;
    return x/(std::exp(x)-1.0);
  }
  
  inline std::shared_ptr<numcxx::DArray1> solve_expfit(int n, double v, double D)
  {
    
    auto h=1.0/(double)(n-1);
    auto pM=numcxx::DSparseMatrix::create(n,n);
    auto pF=numcxx::DArray1::create(n);
    auto pU=numcxx::DArray1::create(n);
    
    auto &M=*pM;
    auto &F=*pF;
    auto &U=*pU;
    
    
    F=0;
    U=0;
    for (int k=0, l=1;k<n-1;k++,l++)
    {
      double g_kl=D* B(v*h/D);
      double g_lk=D* B(-v*h/D);
      M(k,k)+=g_kl/h;
      M(k,l)-=g_kl/h;
      M(l,l)+=g_lk/h;
      M(l,k)-=g_lk/h;
    }

    M(0,0)+=1.0e30;
    M(n-1,n-1)+=1.0e30;
    F(n-1)=1.0e30;

    pM->flush();
    
    auto pUmfpack=numcxx::DSolverUMFPACK::create(pM);
    pUmfpack->update();
    pUmfpack->solve(U,F);
    return pU;
  }


  inline std::shared_ptr<numcxx::DArray1> solve_central(int n, double v, double D)
  {
    
    auto h=1.0/(double)(n-1);
    auto pM=numcxx::DSparseMatrix::create(n,n);
    auto pF=numcxx::DArray1::create(n);
    auto pU=numcxx::DArray1::create(n);
    
    auto &M=*pM;
    auto &F=*pF;
    auto &U=*pU;
    
    
    F=0;
    U=0;
    double g_kl=D - 0.5*(v*h);
    double g_lk=D + 0.5*(v*h);
    for (int k=0, l=1;k<n-1;k++,l++)
    {
      M(k,k)+=g_kl/h;
      M(k,l)-=g_kl/h;
      M(l,l)+=g_lk/h;
      M(l,k)-=g_lk/h;
    }

    M(0,0)+=1.0e30;
    M(n-1,n-1)+=1.0e30;
    F(n-1)=1.0e30;

    pM->flush();
    
    auto pUmfpack=numcxx::DSolverUMFPACK::create(pM);
    pUmfpack->update();
    pUmfpack->solve(U,F);
    return pU;
  }



  inline std::shared_ptr<numcxx::DArray1> solve_upwind(int n, double v, double D)
  {
    
    auto h=1.0/(double)(n-1);
    auto pM=numcxx::DSparseMatrix::create(n,n);
    auto pF=numcxx::DArray1::create(n);
    auto pU=numcxx::DArray1::create(n);
    
    auto &M=*pM;
    auto &F=*pF;
    auto &U=*pU;
    
    
    F=0;
    U=0;
    for (int k=0, l=1;k<n-1;k++,l++)
    {
      double g_kl=D;
      double g_lk=D;
      if (v<0) g_kl-=v*h;
      else  g_lk+=v*h;

      M(k,k)+=g_kl/h;
      M(k,l)-=g_kl/h;
      M(l,l)+=g_lk/h;
      M(l,k)-=g_lk/h;
    }

    M(0,0)+=1.0e30;
    M(n-1,n-1)+=1.0e30;
    F(n-1)=1.0e30;

    pM->flush();
    
    auto pUmfpack=numcxx::DSolverUMFPACK::create(pM);
    pUmfpack->update();
    pUmfpack->solve(U,F);
    return pU;
  }
}


