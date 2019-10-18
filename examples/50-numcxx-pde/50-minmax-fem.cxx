
//(
/// Stationary convergence test for fem method 
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

#include "vtkfigFrame.h"
#include "vtkfigDataSet.h"
#include "vtkfigGridView.h"
#include "vtkfigScalarView.h"
#include "vtkfigXYPlot.h"
#endif

int main(void)
{
  double maxref=9;
  numcxx::Geometry Geometry;
  Geometry.set_points({
        {0,0},
        {0.5,0},
        {1,0},
        {1,1},
        {0.5,1},
        {0,1},
        {0.7,0.7},
    });

  Geometry.set_bfaces({
        {0,1},
        {1,2},
        {2,3},
        {3,4},
        {4,5},
        {5,0}
        });
  
  
  Geometry.set_bfaceregions({1,2,3,4,5,6});
  
  Geometry.set_regionpoints({
      {0.5,0.5}
    });
  Geometry.set_regionnumbers({1});


  auto frame=vtkfig::Frame::New();
  frame->SetSize(800,800);
  frame->SetLayout(2,1);

  auto gridview=vtkfig::GridView::New();
  frame->AddFigure(gridview,0);

  auto solview=vtkfig::ScalarView::New();
  frame->AddFigure(solview,1);
  
  
  double t0;
  double h_intended=0.1;
  double vol= h_intended*h_intended*0.25;
    
  Geometry.set_regionvolumes({vol});
  t0=numcxx::cpu_clock();
  numcxx::SimpleGrid grid(Geometry,"zpaAqDV");
  double t_gen=numcxx::cpu_clock()-t0;

  double hmin,hmax;
  grid.calc_hminmax(hmin,hmax);
  printf("h_intended: %8.3g h_max: %8.3g\n",h_intended,hmax);
  
    
  numcxx::DArray1 bcfac(7);
  numcxx::DArray1 bcval(7);
  bcfac=0;
  bcval=0;
  bcfac(1)=fem2d::Dirichlet;
  bcfac(4)=fem2d::Dirichlet;
  bcval(1)=1;
  bcval(4)=-1;
//  bcfac=fem2d::Dirichlet;;
//  bcval=0;
  
  auto nnodes=grid.npoints();
  auto nodes=grid.get_points();
  
  numcxx::DArray1 source(nnodes);
  numcxx::DArray1 kappa(nnodes);
  kappa=1;
  source=0;
    
  numcxx::DSparseMatrix SGlobal(nnodes,nnodes);
  numcxx::DArray1 Rhs(nnodes);
  numcxx::DArray1 Sol(nnodes);
  
  t0=numcxx::cpu_clock();
  fem2d::assemble_heat_problem(grid,bcfac,bcval,source,kappa,SGlobal, Rhs);
  double t_asm=numcxx::cpu_clock()-t0;
  
  numcxx::DSolverUMFPACK Solver(SGlobal);
  
  t0=numcxx::cpu_clock();
  Solver.update(SGlobal);
  double t_luf=numcxx::cpu_clock()-t0;
  
  t0=numcxx::cpu_clock();
  Solver.solve(Sol,Rhs);

  
#ifdef VTKFIG
  auto griddata=numcxx::vtkfigDataSet(grid);
  griddata->SetPointScalar(Sol ,"Sol");
  
  gridview->SetData(griddata);
  solview->SetData(griddata,"Sol");
  

    
    
  frame->Interact();
#endif

}


