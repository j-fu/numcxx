///
/// \example 46-nonlin-fvm.cxx
///
/// Demo for nonlinear diffusion with fvm
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
      {-2,0},
      {0,0},
      {2,0},
      {2,2},
      {0.5,2},
      {0,1},
      {-0.5,2},
      {-2,2}
    });

  Geometry.set_bfaces({
      {0,1},
      {1,2},
      {2,3},
      {3,4},
      {4,5},
      {5,6},
      {6,7},
      {7,0},
      {5,1},
        });
  
  
  Geometry.set_bfaceregions({1,1,2,3,4,4,3,2,5});
  
  Geometry.set_regionpoints({
      {-0.5,1},
      {0.5,1}
    });
  Geometry.set_regionnumbers({1,2});
  double vol=0.0001;
  Geometry.set_regionvolumes({vol,vol});
  
  
  
  numcxx::SimpleGrid grid(Geometry,"zpaAqDV");
  
  
  
  numcxx::DArray1 bcfac(8);
  numcxx::DArray1 bcval(8);
  bcfac=0;
  bcval=0;
  bcfac(4)=fvm2d::Dirichlet;
  bcfac(1)=fvm2d::Dirichlet;
  bcval(4)=1.0;
  bcval(1)=0.0;

  auto nnodes=grid.npoints();

  numcxx::DArray1 source(nnodes);
  source=0;
  auto fkappa = [](double u, double &kappa,double &dkappa  ) 
    { 
      kappa=0.1+100.0*u*u*u*u;
      dkappa=400.0*u*u*u; 
    };
  

  numcxx::DSparseMatrix SGlobal(nnodes,nnodes);
  numcxx::DSolverUMFPACK Solver(SGlobal);

  numcxx::DArray1 Rhs(nnodes);
  numcxx::DArray1 Sol(nnodes);
  numcxx::DArray1 Res(nnodes);
  numcxx::DArray1 Upd(nnodes);

  Sol=0.0;
  fvm2d::initialize_bc(grid,bcval,Sol);

  

  int iter=0;
  double norm=1.0;
  double oldnorm=1.0;
  double d=1.0;
  double ddelta=1.2;
  
  
  while (iter<100 && norm >1.0e-13)
  {
    fvm2d::assemble_and_apply_nonlinear_heat(grid,bcfac,bcval, source,fkappa, SGlobal,  Sol, Res);
    Solver.update(SGlobal);
    Solver.solve(Upd,Res);
    oldnorm=norm;
    norm=numcxx::norm2(Upd);
    printf("iter=%d norm=%8.4e contract=%8.5e\n", iter,norm,norm/oldnorm);
    Sol=Sol-d*Upd;
    iter++;
    d= std::min(1.0, d*ddelta);
  }


  
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


