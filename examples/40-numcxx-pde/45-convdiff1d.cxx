#include <iostream>
#include <limits>
#include <cassert>
#include <numcxx/numcxx.hxx>
#include <cmath>
#include <iostream>


#include "vtkfigFrame.h"
#include "vtkfigXYPlot.h"




inline double B(double x)
{
  if (std::fabs(x)<1.0e-10) return 1.0;
  return x/(std::exp(x)-1.0);
}
  
void solve_expfit(double v, double D, numcxx::DArray1 &U)
{
  int n=U.size();
  auto h=1.0/(double)(n-1);
  numcxx::DSparseMatrix A(n,n);
  numcxx::DArray1 F(n);
    
  F=0;
  U=0;
  for (int k=0, l=1;k<n-1;k++,l++)
  {
    double g_kl=D* B(v*h/D);
    double g_lk=D* B(-v*h/D);
    A(k,k)+=g_kl/h;
    A(k,l)-=g_kl/h;
    A(l,l)+=g_lk/h;
    A(l,k)-=g_lk/h;
  }

  A(0,0)+=1.0e30;
  A(n-1,n-1)+=1.0e30;
  F(n-1)=1.0e30;
  A.flush();

  numcxx::DSolverUMFPACK Solver(A);
  Solver.update(A);
  Solver.solve(U,F);

}

void solve_central(double v, double D, numcxx::DArray1 &U)
{
  int n=U.size();
  auto h=1.0/(double)(n-1);
  numcxx::DSparseMatrix A(n,n);
  numcxx::DArray1 F(n);
    
  F=0;
  U=0;
  double g_kl=D - 0.5*(v*h);
  double g_lk=D + 0.5*(v*h);
  for (int k=0, l=1;k<n-1;k++,l++)
  {
    A(k,k)+=g_kl/h;
    A(k,l)-=g_kl/h;
    A(l,l)+=g_lk/h;
    A(l,k)-=g_lk/h;
  }
  
  A(0,0)+=1.0e30;
  A(n-1,n-1)+=1.0e30;
  F(n-1)=1.0e30;
  A.flush();
  
  numcxx::DSolverUMFPACK Solver(A);
  Solver.update(A);
  Solver.solve(U,F);

}


void solve_upwind(double v, double D, numcxx::DArray1 &U)
{
  int n=U.size();
  auto h=1.0/(double)(n-1);
  numcxx::DSparseMatrix A(n,n);
  numcxx::DArray1 F(n);
    
  F=0;
  U=0;
  for (int k=0, l=1;k<n-1;k++,l++)
  {
    double g_kl=D;
    double g_lk=D;
    if (v<0) g_kl-=v*h;
    else  g_lk+=v*h;

    A(k,k)+=g_kl/h;
    A(k,l)-=g_kl/h;
    A(l,l)+=g_lk/h;
    A(l,k)-=g_lk/h;
  }

  A(0,0)+=1.0e30;
  A(n-1,n-1)+=1.0e30;
  F(n-1)=1.0e30;
  A.flush();

  numcxx::DSolverUMFPACK Solver(A);
  Solver.update(A);
  Solver.solve(U,F);

}




int main(void)
{

  int n=10;
  
  int maxref=10;
  double D=0.01;
  double v=1.0;
  
  auto frame=vtkfig::Frame::New();
  frame->SetSize(800,800);
  auto plot=vtkfig::XYPlot::New();
  frame->AddFigure(plot);
  
  for (int iref=0;iref<maxref;iref++)
  {
    auto Uexp=numcxx::DArray1(n);
    auto Uupw=numcxx::DArray1(n);
    auto Ucnt=numcxx::DArray1(n);
    auto X=numcxx::linspace(0.0,1.0,n);


    solve_expfit(v,D,Uexp);
    solve_central(v,D,Ucnt);
    solve_upwind(v,D,Uupw);


    frame->SetFrameTitle("n="+std::to_string(n));
    plot->Clear();
    plot->ShowGrid(false);

    plot->SetXRange(0,1);
    plot->SetYRange(-1,1);

    plot->SetLegendPosition(0.3,0.7);
    plot->SetPlotLineType("-");


    plot->SetPlotColor(1,0,0);
    plot->SetPlotLegend("exp");
    plot->AddPlot(*X,Uexp);

    plot->SetPlotColor(0,1,0);
    plot->SetPlotLegend("upw");
    plot->AddPlot(*X,Uupw);

    plot->SetPlotColor(0,0,1);
    plot->SetPlotLegend("cent");
    plot->AddPlot(*X,Ucnt);

    frame->Interact();

    n*=2;
  }

}


