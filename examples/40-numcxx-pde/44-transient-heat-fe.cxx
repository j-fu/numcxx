///
/// \example 40-stationary-heat.cxx
///
/// Heat equation
/// 
#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include <numcxx/numcxx.hxx>
#include <numcxx/simplegrid.hxx>
#include <numcxx/fem2d.hxx>
#ifdef VTKFIG
#include "numcxx/vtkfig-simplegrid.hxx"
#endif

#include "vtkfigFrame.h"
#include "vtkfigDataSet.h"
#include "vtkfigGridView.h"
#include "vtkfigScalarView.h"



int main(void)
{

  numcxx::Geometry Geometry;
  Geometry.set_points({
        {0,0},
        {1,0},
        {1,1},
        {0,1},
        {0.7,0.7},
    });

  Geometry.set_bfaces({
        {0,1},
        {1,2},
        {2,3},
        {3,0}
        });
  
  
  Geometry.set_bfaceregions({1,2,3,4});
  
  Geometry.set_regionpoints({
      {0.5,0.5}
    });

  Geometry.set_regionnumbers({1});
  Geometry.set_regionvolumes({0.0001});
  
  numcxx::SimpleGrid grid(Geometry,"zpaAqDV");


  auto griddata=numcxx::vtkfigDataSet(grid);
  
  auto frame=vtkfig::Frame::New();
  frame->SetSize(800,400);
  frame->SetLayout(2,1);

  auto gridview=vtkfig::GridView::New();
  gridview->SetData(griddata);
  frame->AddFigure(gridview,0);

  auto solview=vtkfig::ScalarView::New();
  solview->SetData(griddata,"Sol");
  solview->SetValueRange(0,1);
  frame->AddFigure(solview,1);
  
  
  numcxx::DArray1 bcfac(6);
  numcxx::DArray1 bcval(6);
  bcfac=0;
  bcval=0;

  bcfac(3)=fem2d::Dirichlet;
  bcfac(1)=fem2d::Dirichlet;
  bcval(3)=1.0;
  bcval(1)=0.0;

  auto nnodes=grid.npoints();

  numcxx::DArray1 source(nnodes);
  numcxx::DArray1 kappa(nnodes);
  kappa=1.0e-2;
  source=0;
  double tau=0.001;
  double theta=0.0;
  double T=100.0;
  bool lump=false;
  int N=T/tau;

  numcxx::DSparseMatrix SGlobal(nnodes,nnodes);
  numcxx::DArray1 Rhs(nnodes);
  numcxx::DArray1 Sol(nnodes);
  numcxx::DArray1 OldSol(nnodes);
  Sol=0.0;
  
  numcxx::DSolverUMFPACK Solver(SGlobal);
 
  for (int n=0;n<N;n++)
  {
    OldSol=Sol;
    fem2d::assemble_transient_heat_problem(grid,bcfac,bcval,source,kappa,tau,theta, lump, OldSol, SGlobal, Rhs);
    
    if (theta>0.0 || !lump)
    {
      Solver.update(SGlobal);
      Solver.solve(Sol,Rhs);
    }
    else
    {
      for (int i=0;i<Sol.size();i++)
       Sol(i)=Rhs(i)/SGlobal(i,i);
    }
    frame->SetFrameTitle("Time="+std::to_string(tau*n));
    
#ifdef VTKFIG
    
    griddata->SetPointScalar(Sol ,"Sol");
    frame->Show();
    
#endif
  }
  frame->Interact();
  
}


