#include <numcxx/numcxx.hxx>
#include <numcxx/simplegrid.hxx>
#include <numcxx/fvm2d.hxx>
#include <cmath>
#include <cassert>
#include <iostream>

namespace fvm2d
{
    
  inline void compute_cell_volume(    
    const int icell,
    const numcxx::DArray2 & points,
    const numcxx::IArray2 & cells,
    double & vol)
  {
    int i0=cells(icell,0);
    int i1=cells(icell,1);
    int i2=cells(icell,2);
    
    // Fill matrix V
    double V00= points(i1,0)- points(i0,0);
    double V01= points(i2,0)- points(i0,0);
    
    double V10= points(i1,1)- points(i0,1);
    double V11= points(i2,1)- points(i0,1);
    
    // Compute determinant
    double det=V00*V11 - V01*V10;
    vol=0.5*det;
  }

  inline void compute_local_formfactors(
    const int icell,
    const numcxx::DArray2 & points,
    const numcxx::IArray2 & cells,
    numcxx::DArray1 & epar,
    numcxx::DArray1 & npar,
    double & vol)
  {
    int i0=cells(icell,0);
    int i1=cells(icell,1);
    int i2=cells(icell,2);
    
    
    // Fill matrix of edge vectors
    double V00= points(i1,0)- points(i0,0);
    double V10= points(i1,1)- points(i0,1);
    
    double V01= points(i2,0)- points(i0,0);
    double V11= points(i2,1)- points(i0,1);
    
    double V02= points(i2,0)- points(i1,0);
    double V12= points(i2,1)- points(i1,1);
    
    
    
    // Compute determinant 
    double det=V00*V11 - V01*V10;
    vol=0.5*det;
    
    double ivol = 1.0/vol;
    
    // squares of edge lengths
    double dd0=V02*V02+V12*V12; // l21
    double dd1=V01*V01+V11*V11; // l20
    double dd2=V00*V00+V10*V10; // l10
    
    
    // contributions to \sigma_kl/h_kl
    epar(0)= (dd1+dd2-dd0)*0.125*ivol;
    epar(1)= (dd2+dd0-dd1)*0.125*ivol;
    epar(2)= (dd0+dd1-dd2)*0.125*ivol;
    
    
    // contributions to \omega_k
    npar(0)= (epar(2)*dd2+epar(1)*dd1)*0.25;
    npar(1)= (epar(0)*dd0+epar(2)*dd2)*0.25;
    npar(2)= (epar(1)*dd1+epar(0)*dd0)*0.25;
  }
  
  
  void assemble_bc(    
    const numcxx::SimpleGrid &grid,
    const numcxx::DArray1& bcfac,
    const numcxx::DArray1& bcval,
    numcxx::DSparseMatrix &SGlobal,
    numcxx::DArray1 &Rhs)
  {
    // Assemble boundary conditions (Dirichlet penalty method)
    int nbfaces=grid.nbfaces();
    auto points=grid.get_points(); // Array of global nodes
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();
    
    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      // Obtain number of boundary condition
      int ireg=bfaceregions(ibface);
      double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
      double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
      double h=sqrt(dx*dx+dy*dy);
      double factor=0.5*h*bcfac(ireg);
      if (bcfac(ireg)>=Dirichlet) factor=bcfac(ireg);
      
      int i;
      i=bfaces(ibface,0);
      SGlobal(i,i)+=factor;
      Rhs(i)+=factor*bcval(ireg);
      
      i=bfaces(ibface,1);
      SGlobal(i,i)+=factor;
      Rhs(i)+=factor*bcval(ireg);
    }

  }
      


  void assemble_apply_bc(    
    const numcxx::SimpleGrid &grid,
    const numcxx::DArray1& bcfac,
    const numcxx::DArray1& bcval,
    numcxx::DSparseMatrix &SGlobal,
    numcxx::DArray1 &Sol,
    numcxx::DArray1 &Rhs)
  {
    // Assemble boundary conditions (Dirichlet penalty method)
    int nbfaces=grid.nbfaces();
    auto points=grid.get_points(); // Array of global nodes
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();
    
    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      // Obtain number of boundary condition
      int ireg=bfaceregions(ibface);
      double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
      double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
      double h=sqrt(dx*dx+dy*dy);
      double factor=0.5*h*bcfac(ireg);
      if (bcfac(ireg)>=Dirichlet) factor=bcfac(ireg);
      
      int i;
      i=bfaces(ibface,0);
      SGlobal(i,i)+=factor;
      Rhs(i)+=factor*(Sol(i)-bcval(ireg));

      
      i=bfaces(ibface,1);
      SGlobal(i,i)+=factor;
      Rhs(i)+=factor*(Sol(i)-bcval(ireg));
    }

  }
      

  void  assemble_heat_problem(
    const numcxx::SimpleGrid &grid,
    const numcxx::DArray1& bcfac,
    const numcxx::DArray1& bcval,
    const numcxx::DArray1& source,
    const numcxx::DArray1& kappa,
    numcxx::DSparseMatrix &SGlobal,
    numcxx::DArray1 &Rhs)
  {
    
    auto ndim=grid.spacedim();      // space dimension
    assert(ndim==2);
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map

    int npoints=grid.npoints();
    int ncells=grid.ncells();
    double vol=0.0;

    // this needs to coorespond to the edge indices
    numcxx::IArray2 edgenodes{{1,2},{0,2},{0,1}};

    const int nnodes_per_cell=3;
    const int nedges_per_cell=3;
    numcxx::DArray1 epar(nnodes_per_cell);
    numcxx::DArray1 npar(nedges_per_cell);
    Rhs=0.0;
    SGlobal(0,0)=0;
    SGlobal.clear();
    
    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {

      compute_local_formfactors(icell, points,cells,epar,npar, vol);
      for (int inode=0;inode<nnodes_per_cell;inode++)
      {
        int k=cells(icell,inode);
        Rhs(k)=source(k)*npar(inode);
      }
      
      for (int iedge=0;iedge<nedges_per_cell;iedge++)
      {
        int k0=cells(icell,edgenodes(iedge,0));
        int k1=cells(icell,edgenodes(iedge,1));
        
        // Assemble fluxes with edge averaged diffusion
        // coefficients into global matrix
        SGlobal(k0,k0)+=epar(iedge)*0.5*(kappa(k0)+kappa(k1));
        SGlobal(k0,k1)-=epar(iedge)*0.5*(kappa(k0)+kappa(k1));
        SGlobal(k1,k0)-=epar(iedge)*0.5*(kappa(k0)+kappa(k1));
        SGlobal(k1,k1)+=epar(iedge)*0.5*(kappa(k0)+kappa(k1));
      }
    }    
    assemble_bc(grid, bcfac, bcval, SGlobal, Rhs);
    SGlobal.flush();
  }
  



  void  assemble_and_apply_nonlinear_heat(
    const numcxx::SimpleGrid &grid,
    const numcxx::DArray1& bcfac,
    const numcxx::DArray1& bcval,
    const numcxx::DArray1& source,
    std::function <void(const double, double&, double&)> fkappa,
    numcxx::DSparseMatrix &SGlobal,
    numcxx::DArray1 &Sol,
    numcxx::DArray1 &Rhs)
  {
    
    auto ndim=grid.spacedim();      // space dimension
    assert(ndim==2);
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map
    
    int npoints=grid.npoints();
    int ncells=grid.ncells();
    double vol=0.0;
    
    // this needs to coorespond to the edge indices
    numcxx::IArray2 edgenodes{{1,2},{0,2},{0,1}};

    const int nnodes_per_cell=3;
    const int nedges_per_cell=3;
    numcxx::DArray1 epar(nnodes_per_cell);
    numcxx::DArray1 npar(nedges_per_cell);
    
    numcxx::DArray1 kappa(nnodes_per_cell);
    numcxx::DArray1 dkappa(nnodes_per_cell);
    

    Rhs=0.0;
    SGlobal(0,0)=0;
    SGlobal.clear();
    

    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {

      compute_local_formfactors(icell, points,cells,epar,npar, vol);
      for (int inode=0;inode<nnodes_per_cell;inode++)
      {
        int k=cells(icell,inode);
        Rhs(k)-=source(k)*npar(inode);
        fkappa(Sol(k),kappa(inode),dkappa(inode));
      }

      for (int iedge=0;iedge<nedges_per_cell;iedge++)
      {
        int i0=edgenodes(iedge,0);
        int i1=edgenodes(iedge,1);
        int k0=cells(icell,i0);
        int k1=cells(icell,i1);
        
        
        double flux=epar(iedge)*0.5*(kappa(i0)+kappa(i1))*(Sol(k0)-Sol(k1));
        Rhs(k0)+=flux;
        Rhs(k1)-=flux;
        
        SGlobal(k0,k0)+= epar(iedge)*( 0.5*(kappa(i0)+kappa(i1))+0.5*dkappa(i0)*(Sol(k0)-Sol(k1)));
        SGlobal(k0,k1)+= epar(iedge)*(-0.5*(kappa(i0)+kappa(i1))+0.5*dkappa(i1)*(Sol(k0)-Sol(k1)));
        SGlobal(k1,k0)+= epar(iedge)*(-0.5*(kappa(i0)+kappa(i1))-0.5*dkappa(i0)*(Sol(k0)-Sol(k1)));
        SGlobal(k1,k1)+= epar(iedge)*( 0.5*(kappa(i0)+kappa(i1))-0.5*dkappa(i1)*(Sol(k0)-Sol(k1)));
      }
    }    
    assemble_apply_bc(grid, bcfac, bcval, SGlobal, Sol,Rhs);
    SGlobal.flush();
  }
  

  void initialize_bc(
    numcxx::SimpleGrid &grid,// Discretization grid
    numcxx::DArray1& g,
    numcxx::DArray1& Sol
    )
  {
    
    int nbfaces=grid.nbfaces();
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();
    
    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      int i;
      int ireg=bfaceregions(ibface);
      
      i=bfaces(ibface,0);
      Sol(i)=g(ireg);
      
      i=bfaces(ibface,1);
      Sol(i)=g(ireg);
    }
    
  }



  double l2norm(const numcxx::SimpleGrid &grid, 
                const numcxx::DArray1 &u)
  {
    auto ndim=grid.spacedim();
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map
    int npoints=grid.npoints();
    int ncells=grid.ncells();

    double vol=0.0;


    numcxx::DArray1 epar(ndim+1);
    numcxx::DArray1 npar(ndim+1);

    
    double norm=0.0;
    for (int icell=0; icell<ncells; icell++)
    {
      double vol;
      compute_local_formfactors(icell, points,cells,epar,npar, vol);

      for (int i=0;i<=ndim;i++)
        norm+=u(cells(icell,i))*u(cells(icell,i))*npar(i);
    }
    return sqrt(norm);
  }

  double h1norm(const numcxx::SimpleGrid &grid, 
                const numcxx::DArray1 &u)
  {
    auto ndim=grid.spacedim();
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map
    int npoints=grid.npoints();
    int ncells=grid.ncells();

    double vol=0.0;


    numcxx::DArray1 epar(ndim+1);
    numcxx::DArray1 npar(ndim+1);

    
    double norm=0.0;
    for (int icell=0; icell<ncells; icell++)
    {
      double vol;
      compute_local_formfactors(icell, points,cells,epar,npar, vol);
      
      double d01=u(cells(icell,0))-u(cells(icell,1));
      double d02=u(cells(icell,0))-u(cells(icell,2));
      double d12=u(cells(icell,1))-u(cells(icell,2));


      norm+=epar(0)*d12*d12+epar(1)*d02*d02+epar(2)*d01*d01;
    }


    return sqrt(norm);
  }






////////////////////////////////////////////////////////

  





}
