///
/// Multigrid test, under development
/// 
#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include <numcxx/numcxx.hxx>
#include <numcxx/simplegrid.hxx>
#include <numcxx/fvm2d.hxx>
#include <numcxx/tarray.hxx>
#include <netlib/netlib.hxx>

#ifdef VTKFIG
#include "numcxx/vtkfig-simplegrid.hxx"
#endif

#include "vtkfigFrame.h"
#include "vtkfigDataSet.h"
#include "vtkfigGridView.h"
#include "vtkfigScalarView.h"


template<typename T> 
class TPreconMG: public  numcxx::TLinSolver<T>
{
public:
  int nlev;
  std::vector<std::shared_ptr<numcxx::SimpleGrid>> Grids;
  std::vector<std::shared_ptr<numcxx::TPreconJacobi<T>>> Smoothers;
  std::vector<std::shared_ptr<numcxx::TSparseMatrix<T>>> Matrices;
  std::shared_ptr<numcxx::TSolverUMFPACK<T>> CoarseSolver;
  std::vector<int> ndir; // number of points in coordinate directions
  std::vector<int> nnodes; // numbers of points on grid
  std::vector<std::shared_ptr<numcxx::TArray1<T> >> Sol;
  int nsmooth=1; 
  int gamma=1; // coarse grid corr.
  double cdamp=1;
  double sdamp=1;
  
  
  void prolongate(const int lev, const numcxx::DArray1 &ucoarse,  numcxx::DArray1 &ufine) const
  {
    assert(ucoarse.shape(0)==nnodes[lev-1]);
    assert(ufine.shape(0)==nnodes[lev]);

    int nyc=ndir[lev-1];
    int nxc=ndir[lev-1];
    int nxf=ndir[lev];
    int nyf=ndir[lev];
    
    ufine=0.0;
    int ifp=0;

    for (int iyc=0;iyc<nyc-1;iyc++, ifp+=2*nxc)
      for (int ixc=0;ixc<nxc-1;ixc++, ifp+=2)
      {
        
	int icp=ixc+iyc*nxc;
        // Interpolate data from coarse grid rectangle
        // to four fine grid rectangles

        // Corners
        ufine(ifp)=ucoarse(icp);
        ufine(ifp+2)=ucoarse(icp+1);
        ufine(ifp+2*nxf)=ucoarse(icp+nxc);
        ufine(ifp+2*nxf+2)=ucoarse(icp+nxc+1);

        
        // Edges
        ufine(ifp+nxf)=0.5*(ucoarse(icp)+ucoarse(icp+nxc));
        ufine(ifp+nxf+2)=0.5*(ucoarse(icp+1)+ucoarse(icp+nxc+1));
        ufine(ifp+1)=0.5*(ucoarse(icp)+ucoarse(icp+1));
        ufine(ifp+2*nxf+1)=0.5*(ucoarse(icp+nxc)+ucoarse(icp+nxc+1));

        // Midpoint
        ufine(ifp+nxf+1)=0.25*(ucoarse(icp)+ucoarse(icp+nxc)+ucoarse(icp+1)+ucoarse(icp+nxc+1));
      }

  }

  void restrict(const int lev, numcxx::DArray1 &ucoarse,  const numcxx::DArray1 &ufine) const
  {
    assert(ucoarse.shape(0)==nnodes[lev-1]);
    assert(ufine.shape(0)==nnodes[lev]);
    
    int nyc=ndir[lev-1];
    int nxc=ndir[lev-1];
    int nxf=ndir[lev];
    int nyf=ndir[lev];
    
    ucoarse=0.0;
    int ifp=0;
    // Corner loop
    for (int iyc=0;iyc<nyc-1;iyc++, ifp+=2*nxc)
      for (int ixc=0;ixc<nxc-1;ixc++, ifp+=2)
      {
        
	int icp=ixc+iyc*nxc;
        ucoarse(icp)=ufine(ifp);
        ucoarse(icp+1)=ufine(ifp+2);
        ucoarse(icp+nxc)=ufine(ifp+2*nxf);
        ucoarse(icp+nxc+1)=ufine(ifp+2*nxf+2);
      }


    ifp=0;
    for (int iyc=0;iyc<nyc-1;iyc++, ifp+=2*nxc)
      for (int ixc=0;ixc<nxc-1;ixc++, ifp+=2)
      {

	int icp=ixc+iyc*nxc;
        double fac;
        // Restriction from edges
        // south edge
        ucoarse(icp)+=0.5*ufine(ifp+1);
        ucoarse(icp+1)+=0.5*ufine(ifp+1);

        // west edge
        ucoarse(icp)+=0.5*ufine(ifp+nxf);
        ucoarse(icp+nxc)+=0.5*ufine(ifp+nxf);

        // east edge  - we need to do this only once
        if (ixc==(nxc-2))
        {
          ucoarse(icp+1)+=0.5*ufine(ifp+nxf+2);
          ucoarse(icp+nxc+1)+=0.5*ufine(ifp+nxf+2);
        }

        // north edge - we need to do this only once
        if (iyc==(nyc-2))
        {
          ucoarse(icp+nxc)+=0.5*ufine(ifp+2*nxf+1);
          ucoarse(icp+nxc+1)+=0.5*ufine(ifp+2*nxf+1);
        }

        // Restriction from midpoint
        ucoarse(icp)+=0.25*ufine(ifp+nxf+1);
        ucoarse(icp+1)+=0.25*ufine(ifp+nxf+1);
        ucoarse(icp+nxc)+=0.25*ufine(ifp+nxf+1);
        ucoarse(icp+nxc+1)+=0.25*ufine(ifp+nxf+1);
      }
  }

  void solve(numcxx::TArray1<T> & Sol, const numcxx::TArray1<T> & Rhs) const
  {
    Sol.resize(Rhs.size());
    Sol=0.0;
    cycle(nlev,Sol, Rhs);
  }

  void cycle(int ilev, numcxx::TArray1<T> & Sol,  const numcxx::TArray1<T> & Rhs)  const
  {
    if (ilev==0)
    {
      CoarseSolver->solve(Sol,Rhs);
    }
    else
    {
      auto Res=numcxx::DArray1(nnodes[ilev]);
      auto Upd=numcxx::DArray1(nnodes[ilev]);
      auto & A=*(Matrices[ilev]);
      auto & M=*(Smoothers[ilev]);
      for (int i=0;i<nsmooth;i++)
      {
        Res=A*Sol-Rhs;
        M.solve(Upd,Res);
        Sol=Sol-sdamp*Upd;
      }

      Res=A*Sol-Rhs;
      auto CSol=numcxx::DArray1(nnodes[ilev-1]);
      auto CRhs=numcxx::DArray1(nnodes[ilev-1]);
      CSol=0.0;
      restrict(ilev,CRhs,Res);
      for (int i=0; i<gamma; i++)
        cycle(ilev-1, CSol, CRhs);
      prolongate(ilev,CSol, Upd);
      Sol=Sol-cdamp*Upd;
      for (int i=0;i<nsmooth; i++)
      {
        Res=A*Sol-Rhs;
        M.solve(Upd,Res);
        Sol=Sol-sdamp*Upd;
      }
    }
  }
  
};
  

int main(void)
{
  TPreconMG<double> mg;
  int maxlev=2;
  int n0=10;
  mg.nsmooth=1;
  mg.cdamp=1;
  mg.sdamp=1;
  mg.gamma=1;
  int max_iter=40;
  bool CG=false;
  
  auto frame=vtkfig::Frame::New();
  frame->SetSize(800,400);
  frame->SetLayout(2,1);

  auto gridview=vtkfig::GridView::New();
  frame->AddFigure(gridview,0);
  
  auto solview=vtkfig::ScalarView::New();
  frame->AddFigure(solview,1);


  std::shared_ptr<numcxx::SimpleGrid> grid;
  std::shared_ptr<numcxx::DArray1> Rhs;
  std::shared_ptr<numcxx::DArray1> Exact;
  std::shared_ptr<numcxx::DSparseMatrix> SGlobal;
  int nnodes;
  
  for (int ilev=0; ilev<maxlev;ilev++, n0*=2)
  {
    auto X=numcxx::linspace(0.0,1.0,n0+1);
    auto Y=numcxx::linspace(0.0,1.0,n0+1);
    
    grid=std::make_shared<numcxx::SimpleGrid>(*X,*Y);
    mg.Grids.push_back(grid);
    mg.ndir.push_back(n0+1);
    std::cout << n0 << " " <<  grid->npoints() << std::endl;
    auto sol=numcxx::DArray1::create(grid->npoints());
    mg.Sol.push_back(sol);
    
    
    numcxx::DArray1 bcfac(6);
    numcxx::DArray1 bcval(6);
    bcfac=0;
    bcval=0;
    bcfac(3)=fvm2d::Dirichlet;
    bcfac(1)=fvm2d::Dirichlet;
    bcfac(2)=fvm2d::Dirichlet;
    bcfac(4)=fvm2d::Dirichlet;

    nnodes=grid->npoints();
    mg.nnodes.push_back(nnodes);
    
    numcxx::DArray1 source(nnodes);
    numcxx::DArray1 kappa(nnodes);
    numcxx::DArray1 exact(nnodes);
    Exact=numcxx::DArray1::create(nnodes);

    kappa=1;
    source=0;
    auto nodes=grid->get_points();

    for (int i=0; i<nnodes;i++)
    {
      double x=nodes(i,0);
      double y=nodes(i,1);
      (*Exact)[i]=sin(M_PI*x)*sin(M_PI*y);
      source[i]=2.0*M_PI*M_PI*(*Exact)[i];
    }

    
    SGlobal=numcxx::DSparseMatrix::create(nnodes,nnodes);
    Rhs=numcxx::DArray1::create(nnodes);
  
    fvm2d::assemble_heat_problem(*grid,bcfac,bcval,source,kappa,*SGlobal, *Rhs);

    mg.Matrices.push_back(SGlobal);

    if (ilev==0)
    {
      mg.CoarseSolver=numcxx::DSolverUMFPACK::create(SGlobal);
      mg.CoarseSolver->update(*SGlobal);
    }
    else if (0)
    {
      auto griddata=numcxx::vtkfigDataSet(*(mg.Grids)[ilev-1]);
      *(mg.Sol)[ilev-1]=1.0;
      mg.prolongate(ilev, *(mg.Sol)[ilev-1],*(mg.Sol)[ilev]);
      mg.restrict(ilev, *(mg.Sol)[ilev-1],*(mg.Sol)[ilev]);
      griddata->SetPointScalar(*(mg.Sol)[ilev-1] ,"Sol");
      gridview->SetData(griddata);
      solview->SetData(griddata,"Sol");
      frame->Interact();
    
    }
    mg.Smoothers.push_back(numcxx::DPreconJacobi::create(SGlobal));
    mg.Smoothers[ilev]->update(*SGlobal);
    mg.nlev=ilev;
  }

  auto Solver=numcxx::DSolverUMFPACK::create(SGlobal);
  Solver->solve(*Exact,*Rhs);

  
  numcxx::DArray1 Sol(nnodes);
  numcxx::DArray1 Res(nnodes);
  numcxx::DArray1 Upd(nnodes);
  Sol=0.0;
  Upd=0.0;
  double oldnorm=0.0;
  if (not CG)
  {
    for (int nsteps=0;nsteps<max_iter;nsteps++)
    {
      Res=(*SGlobal)*Sol-(*Rhs);
      mg.solve(Upd,Res);
      Sol=Sol-Upd;
      double norm=numcxx::norm2(Sol-(*Exact));
      if (nsteps>0)
        printf("norm=%g contr=%g\n", norm,norm/oldnorm);
      else
        printf("norm=%g\n", norm);
      oldnorm=norm;
    }
  }
  else
  {
    double eps=1.0e-10;
    int retval=netlib::CG(*SGlobal,Sol, *Rhs,mg, max_iter,eps);
    double norm=numcxx::norm2(Sol-(*Exact));
    printf("rv= %d max_iter=%d eps=%g norm=%g\n",retval, max_iter, eps,norm);
  }
  
  
  
  auto griddata=numcxx::vtkfigDataSet(*grid);
  griddata->SetPointScalar(Sol ,"Sol");
  gridview->SetData(griddata);
  solview->SetData(griddata,"Sol");
  frame->Interact();
  
  
}


