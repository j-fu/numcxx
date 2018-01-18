///
/// \example 42-convtest-fem-sin.cxx
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


  auto frame=vtkfig::Frame::New();
  frame->SetSize(800,800);
  frame->SetLayout(2,2);

  auto gridview=vtkfig::GridView::New();
  frame->AddFigure(gridview,0);

  auto solview=vtkfig::ScalarView::New();
  frame->AddFigure(solview,1);
  
  auto errplot=vtkfig::XYPlot::New();
  frame->AddFigure(errplot,2);
  errplot->SetXTitle("log10(h)");
  errplot->SetYTitle("log10(error)");

  auto timeplot=vtkfig::XYPlot::New();
  frame->AddFigure(timeplot,3);
  timeplot->SetXTitle("log10(N)");
  timeplot->SetYTitle("log10(time)");




  std::vector<double> H;
  std::vector<double> N;
  std::vector<double> TAsm;
  std::vector<double> TLuf;
  std::vector<double> TLus;
  std::vector<double> TGen;;


  std::vector<double> L2Error;
  std::vector<double> H1Error;
  
  
  for (int iref=0;iref<maxref;iref++)
  {
    double t0;
    double h_intended=0.5*pow(2.0,-iref);
    double vol= h_intended*h_intended*0.25;
    
    Geometry.set_regionvolumes({vol});
    t0=numcxx::cpu_clock();
    numcxx::SimpleGrid grid(Geometry,"zpaAqDV");
    double t_gen=numcxx::cpu_clock()-t0;

    double hmin,hmax;
    grid.calc_hminmax(hmin,hmax);
    printf("h_intended: %8.3g h_max: %8.3g\n",h_intended,hmax);
    
    
    numcxx::DArray1 bcfac(6);
    numcxx::DArray1 bcval(6);
    bcfac=0;
    bcval=0;
    
    bcfac(1)=fem2d::Dirichlet;
    bcfac(2)=fem2d::Dirichlet;
    bcfac(3)=fem2d::Dirichlet;
    bcfac(4)=fem2d::Dirichlet;
    bcval(1)=0.0;
    bcval(2)=0.0;
    bcval(3)=0.0;
    bcval(4)=0.0;
    
    auto nnodes=grid.npoints();
    auto nodes=grid.get_points();
    
    numcxx::DArray1 source(nnodes);
    numcxx::DArray1 exact(nnodes);
    numcxx::DArray1 kappa(nnodes);
    kappa=1;
    
    for (int i=0; i<nnodes;i++)
    {
      double x=nodes(i,0);
      double y=nodes(i,1);
      exact[i]=sin(M_PI*x)*sin(M_PI*y);
      source[i]=2.0*M_PI*M_PI*exact[i];
    }
    
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
    double t_lus=numcxx::cpu_clock()-t0;
    
    double l2error=fem2d::l2norm(grid,Sol-exact);
    double h1error=fem2d::h1norm(grid,Sol-exact);

    printf("h: %8.3e l2: %8.3e h1: %8.3e\n",hmax,l2error,h1error);
    H.push_back(log10(hmax));
    H1Error.push_back(log10(h1error));
    L2Error.push_back(log10(l2error));

    printf("time/ms gen: %8.3f asm: %8.3f luf: %8.3f lus: %8.3f\n",
           t_gen*1000.0, t_asm*1000.0, t_luf*1000.0, t_lus*1000.0);
    
    N.push_back(log10(grid.npoints()));;
    TAsm.push_back(log10(t_asm));
    TLuf.push_back(log10(t_luf));
    TLus.push_back(log10(t_lus));
    TGen.push_back(log10(t_gen));


#ifdef VTKFIG
    auto griddata=numcxx::vtkfigDataSet(grid);
    griddata->SetPointScalar(Sol ,"Sol");

    gridview->SetData(griddata);
    solview->SetData(griddata,"Sol");

    frame->SetFrameTitle("iref="+std::to_string(iref));
    if (iref>0)
    {
      double x0=-3; 
      double x1=0;
      auto X=numcxx::DArray1{x0,x1};
      auto OH1=numcxx::DArray1{x0,x1};
      auto OH2=numcxx::DArray1{2*x0,2*x1};
      
      errplot->Clear();


      errplot->SetLegendPosition(0.7,0.4);
      errplot->SetXRange(x0,x1);
      errplot->SetYRange(-6,0);


      errplot->SetPlotLineType("-");
      errplot->SetPlotMarkerType("o");
      errplot->SetMarkerSize(2);


      errplot->SetPlotColor(0.2,0.5,0.2);
      errplot->SetPlotLegend("O(h**2)");
      errplot->AddPlot(X,OH2);
      errplot->SetPlotColor(0.5,0.2,0.2);
      errplot->SetPlotLegend("O(h)");
      errplot->AddPlot(X,OH1);

      errplot->SetPlotColor(0.0,0.5,0);
      errplot->SetPlotLegend("L2");

      errplot->AddPlot(H,L2Error);
      errplot->SetPlotLegend("H1");
      errplot->SetPlotColor(0.5,0,0);
      errplot->AddPlot(H,H1Error);

      
      double n0=0;
      double n1=7;
      double n11=6;
      double t0=-5;
      double t1=2;

      auto XN=numcxx::DArray1{n0,n1};
      auto XN1=numcxx::DArray1{n0,n11};
      auto ON1=numcxx::DArray1{n0-6.0,n1-6.0};
      auto ON32=numcxx::DArray1{1.5*n0-8.0, 1.5*n11-8.0};


      timeplot->SetLegendPosition(0.3,0.6);
      timeplot->SetLegendSize(0.15,0.3);
      timeplot->SetXRange(n0,n1);
      timeplot->SetYRange(t0,t1);


      timeplot->Clear();
      timeplot->SetPlotLineType("-");
      timeplot->SetPlotMarkerType("o");
      timeplot->SetMarkerSize(2);


      timeplot->SetPlotLegend("Gen");
      timeplot->SetPlotColor(0.5,0,0);
      timeplot->AddPlot(N,TGen);

      timeplot->SetPlotLegend("Asm");
      timeplot->SetPlotColor(0.0,0.5,0);
      timeplot->AddPlot(N,TAsm);

      timeplot->SetPlotLegend("LU fact");
      timeplot->SetPlotColor(0.0,0.0,0.5);
      timeplot->AddPlot(N,TLuf);

      timeplot->SetPlotLegend("LU solve");
      timeplot->SetPlotColor(0.25,0.25,0);
      timeplot->AddPlot(N,TLus);

      timeplot->SetPlotLegend("O(N**(3/2))");
      timeplot->SetPlotColor(0.0,0.0,1.0);
      timeplot->AddPlot(XN1,ON32);

      timeplot->SetPlotLegend("O(N)");
      timeplot->SetPlotColor(1.0,0,0);
      timeplot->AddPlot(XN,ON1);;


    }
    frame->Interact();
#endif
  }

}


