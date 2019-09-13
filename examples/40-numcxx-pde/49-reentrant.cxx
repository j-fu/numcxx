///
/// \example 41-stationary-heat-fv.cxx
///
/// Finite volume method for stationary heat equation
/// 
#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include <numcxx/numcxx.hxx>
#include <numcxx/simplegrid.hxx>
#include <numcxx/fvm2d.hxx>
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
      {0.0,0.0},
      {0.5,0.0},
      {1.0,0.0},
      {1.0,0.5},
      {0.5,0.5},
      {0.5,1.0},
      {0.0,1.0},
    });

  Geometry.set_bfaces({
      {0,1},
      {1,2},
      {2,3},
      {3,4},
      {4,5},
      {5,6},
      {6,0},
        });
  
  
  Geometry.set_bfaceregions({1,2,3,4,5,6,7});
  
  Geometry.set_regionpoints({
      {0.1, 0.1}
    });
  Geometry.set_regionnumbers({1});
  Geometry.set_regionvolumes({0.001});
  
  
  
  numcxx::SimpleGrid grid(Geometry,"zpaAqDV");
  
  
  
  numcxx::DArray1 bcfac(9);
  numcxx::DArray1 bcval(9);
  bcfac=0;
  bcval=0;
  bcfac(4)=fvm2d::Dirichlet;
  bcval(4)=1.0;

  bcfac(5)=fvm2d::Dirichlet;
  bcval(5)=1.0;

  bcfac(1)=fvm2d::Dirichlet;
  bcval(1)=0.0;

  bcfac(7)=fvm2d::Dirichlet;
  bcval(7)=0.0;


  auto nnodes=grid.npoints();

  numcxx::DArray1 source(nnodes);
  numcxx::DArray1 kappa(nnodes);
  kappa=1;
  source=0;

  numcxx::DSparseMatrix SGlobal(nnodes,nnodes);
  numcxx::DArray1 Rhs(nnodes);
  numcxx::DArray1 Sol(nnodes);
  
  fvm2d::assemble_heat_problem(grid,bcfac,bcval,source,kappa,SGlobal, Rhs);
  numcxx::DSolverUMFPACK Solver(SGlobal);
  
  Solver.update(SGlobal);
  Solver.solve(Sol,Rhs);
  
  
#ifdef VTKFIG
  auto griddata=numcxx::vtkfigDataSet(grid);
  griddata->SetPointScalar(Sol ,"Sol");

  auto frame=vtkfig::Frame::New();
  frame->SetSize(800,400);
  frame->SetLayout(2,1);

  auto gridview=vtkfig::GridView::New();
  gridview->SetData(griddata);
  frame->AddFigure(gridview,0);
  
  auto solview=vtkfig::ScalarView::New();
  solview->SetData(griddata,"Sol");
  frame->AddFigure(solview,1);


  
  frame->Interact();
#endif
}


