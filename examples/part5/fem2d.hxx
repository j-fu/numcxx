#ifndef FEM2D_H
#define FEM2D_H

#include <numcxx/simplegrid.hxx>
#include <numcxx/tsparsematrix.hxx>
#include <cmath>
#include <iostream>

namespace fem2d
{

  const double Dirichlet=1.0e30;

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

  
  inline void compute_local_stiffness_matrix(
    const int icell,
    const numcxx::DArray2 & points,
    const numcxx::IArray2 & cells,
    numcxx::DArray2 & SLocal,
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
    double fac = 0.5/det;
    
    // Compute entries of local stiffness matrix
    SLocal(0,0)= fac * (  ( V10-V11 )*( V10-V11 )+( V01-V00 )*( V01-V00 ) );
    SLocal(0,1)= fac * (  ( V10-V11 )* V11          - ( V01-V00 )*V01 );
    SLocal(0,2)= fac * ( -( V10-V11 )* V10          + ( V01-V00 )*V00 );
    
    SLocal(1,1)=  fac * (  V11*V11 + V01*V01 );
    SLocal(1,2)=  fac * ( -V11*V10 - V01*V00 );
    
    SLocal(2,2)=  fac * ( V10*V10+ V00*V00 );
    
    SLocal(1,0)=SLocal(0,1);
    SLocal(2,0)=SLocal(0,2);
    SLocal(2,1)=SLocal(1,2);

    vol=0.5*det;
  }





  inline void
  assemble_simple_heat_problem(
    std::shared_ptr<numcxx::SimpleGrid> &pGrid,
    std::shared_ptr<numcxx::DArray1> pBCfac,
    std::shared_ptr<numcxx::DArray1> pBCval,
    std::shared_ptr<numcxx::DSparseMatrix>&pS,
    std::shared_ptr<numcxx::DArray1> &pRhs
    )
    {
    auto & grid=*pGrid;
    auto & bcfac=*pBCfac;
    auto & bcval=*pBCval;
    auto & Rhs=*pRhs;
    auto & S=*pS;
    
    auto ndim=grid.spacedim();
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map

    int npoints=grid.npoints();
    int ncells=grid.ncells();
    double vol=0.0;

    // Local stiffness matrix
    auto pSLocal=numcxx::DArray2::create(ndim+1, ndim+1);
    auto &SLocal=*pSLocal;
    
    Rhs=0.0;
    S(0,0)=0;
    S.clear();

    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {
      compute_local_stiffness_matrix(icell, points,cells, SLocal, vol);
      // Assemble into global stiffness matrix
      for (int i=0;i<=ndim;i++)
        for (int j=0;j<=ndim;j++)
          S(cells(icell,i),cells(icell,j))+=SLocal(i,j);
    }    



    // Assemble boundary conditions (Dirichlet penalty method)
    int nbfaces=grid.nbfaces();
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();
    
    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      // Obtain number of boundary condition
      int ireg=bfaceregions(ibface);
      int i0=bfaces(ibface,0);
      int i1=bfaces(ibface,1);
      
      // Check if it is "Dirichlet"
      if (bcfac(ireg)>=Dirichlet)
      {
        // Assemble penalty values
        S(i0,i0)+=bcfac(ireg);
        Rhs(i0)+=Dirichlet*bcval(ireg);
        
        S(i1,i1)+=bcfac(ireg);
        Rhs(i1)+=Dirichlet*bcval(ireg);
        
      }
      else if (bcfac(ireg)>0.0)
      {

        double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
        double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
        double h=sqrt(dx*dx+dy*dy);
        
        S(i0,i0)+=h*bcfac(ireg)/3.0;
        S(i0,i1)+=h*bcfac(ireg)/6.0;
        
        S(i1,i0)+=h*bcfac(ireg)/6.0;
        S(i1,i1)+=h*bcfac(ireg)/3.0;
        
        Rhs(i0)+=0.5*h*bcval(ireg);
        Rhs(i0)+=0.5*h*bcval(ireg);
      }
    }
    S.flush();
    }




  inline void
  assemble_heat_problem_with_source(
    std::shared_ptr<numcxx::SimpleGrid> &pGrid,
    std::shared_ptr<numcxx::DArray1> pBCfac,
    std::shared_ptr<numcxx::DArray1> pBCval,
    std::shared_ptr<numcxx::DArray1> pSource,
    std::shared_ptr<numcxx::DArray1> pKappa,
    std::shared_ptr<numcxx::DSparseMatrix>&pS,
    std::shared_ptr<numcxx::DArray1> &pRhs
    )
    {
    auto & grid=*pGrid;
    auto & bcfac=*pBCfac;
    auto & bcval=*pBCval;
    auto & Rhs=*pRhs;
    auto & S=*pS;
    auto & source=*pSource;
    auto & kappa=*pKappa;

    auto ndim=grid.spacedim();
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map

    int npoints=grid.npoints();
    int ncells=grid.ncells();
    double vol=0.0;

    // Local stiffness matrix
    auto pSLocal=numcxx::DArray2::create(ndim+1, ndim+1);
    auto &SLocal=*pSLocal;
    auto pMLocal0=numcxx::DMatrix::create({{2,1,1},{1,2,1},{1,1,2}}); // from Stroud quadrature...
    auto &MLocal0=*pMLocal0;
    MLocal0*=1.0/12.0;
    
    Rhs=0.0;
    S(0,0)=0;
    S.clear();
    double third=1.0/3.0;

    // Loop over all elements (cells) of the triangulation
    for (int icell=0; icell<ncells; icell++)
    {
      compute_local_stiffness_matrix(icell, points,cells, SLocal, vol);
      // Assemble into global stiffness matrix
      double klocal=(kappa(cells(icell,0))+kappa(cells(icell,1))+kappa(cells(icell,2)))/3.0;
      for (int i=0;i<=ndim;i++)
        for (int j=0;j<=ndim;j++)
        {
          Rhs(cells(icell,i))+=vol*MLocal0(i,j)*source(cells(icell,j));
          S(cells(icell,i),cells(icell,j))+=klocal*SLocal(i,j);
        }
    }    
    


    // Assemble boundary conditions (Dirichlet penalty method)
    int nbfaces=grid.nbfaces();
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();
    
    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      // Obtain number of boundary condition
      int ireg=bfaceregions(ibface);
      int i0=bfaces(ibface,0);
      int i1=bfaces(ibface,1);
      
      // Check if it is "Dirichlet"
      if (bcfac(ireg)>=Dirichlet)
      {
        // Assemble penalty values
        S(i0,i0)+=bcfac(ireg);
        Rhs(i0)+=bcfac(ireg)*bcval(ireg);
        
        S(i1,i1)+=bcfac(ireg);
        Rhs(i1)+=bcfac(ireg)*bcval(ireg);
        
      }
      else if (bcfac(ireg)>0.0)
      {

        double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
        double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
        double h=sqrt(dx*dx+dy*dy);
        
        S(i0,i0)+=h*bcfac(ireg)/3.0;
        S(i0,i1)+=h*bcfac(ireg)/6.0;
        
        S(i1,i0)+=h*bcfac(ireg)/6.0;
        S(i1,i1)+=h*bcfac(ireg)/3.0;
        
        Rhs(i0)+=0.5*h*bcval(ireg);
        Rhs(i1)+=0.5*h*bcval(ireg);
      }
    }
    S.flush();
    }
  


  
/////////////////////////////////////////////////////////
  void assemble_transient_heat_matrix_and_rhs(
    std::shared_ptr<numcxx::SimpleGrid> &pGrid,
    std::shared_ptr<numcxx::DArray1> pBCfac,
    std::shared_ptr<numcxx::DArray1> pBCval,
    std::shared_ptr<numcxx::DArray1> pSource,
    std::shared_ptr<numcxx::DArray1> pKappa,
    double tau, // time step
    double theta, //  choice of method
    bool lump, //  lump mass matrix
    std::shared_ptr<numcxx::DArray1> pOldSol,
    std::shared_ptr<numcxx::DSparseMatrix>&pS,
    std::shared_ptr<numcxx::DArray1> &pRhs
    )
  {
    auto & grid=*pGrid;
    auto & bcfac=*pBCfac;
    auto & bcval=*pBCval;
    auto & Rhs=*pRhs;
    auto & OldSol=*pOldSol;
    auto & S=*pS;
    auto & source=*pSource;
    auto & kappa=*pKappa;

    auto ndim=grid.spacedim();
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map
    int npoints=grid.npoints();
    int ncells=grid.ncells();
    double tauinv=1.0/tau;
    
    // Local stiffness matrix
    auto pSLocal=numcxx::DMatrix::create(ndim+1, ndim+1);
    auto &SLocal=*pSLocal;
    
    // Local mass matrix
   
    auto pMFull=numcxx::DMatrix::create({{2,1,1},{1,2,1},{1,1,2}}); // from Stroud quadrature...
    auto pMLumped=numcxx::DMatrix::create({{4,0,0},{0,4,0},{0,0,4}}); // mass lumping
    std::shared_ptr<numcxx::DMatrix> pMLocal0;
    if (lump)
      pMLocal0=pMLumped;
    else
      pMLocal0=pMFull;
    
    auto &MLocal0=*pMLocal0;
    MLocal0*=1.0/12.0;
    
    // Local matrix of coordinate differences
    auto pV=numcxx::DMatrix::create(ndim, ndim);
    auto &V=*pV;
    
    // Loop over all elements (cells) of the triangulation
    Rhs=0.0;
    S(0,0)=0;
    S.clear();
    for (int icell=0; icell<ncells; icell++)
    {
      // Fill matrix V
      V(0,0)= points(cells(icell,1),0)- points(cells(icell,0),0);
      V(0,1)= points(cells(icell,2),0)- points(cells(icell,0),0);
      
      V(1,0)= points(cells(icell,1),1)- points(cells(icell,0),1);
      V(1,1)= points(cells(icell,2),1)- points(cells(icell,0),1);
      

      // Compute determinant
      double det=V(0,0)*V(1,1) - V(0,1)*V(1,0);
      double vol = 0.5*det;
      double invvol=1.0/vol;
      
      // Compute entries of local stiffness matrix
      SLocal(0,0)= invvol * (  ( V(1,0)-V(1,1) )*( V(1,0)-V(1,1) )+( V(0,1)-V(0,0) )*( V(0,1)-V(0,0) ) );
      SLocal(0,1)= invvol * (  ( V(1,0)-V(1,1) )* V(1,1)          - ( V(0,1)-V(0,0) )*V(0,1) );
      SLocal(0,2)= invvol * ( -( V(1,0)-V(1,1) )* V(1,0)          + ( V(0,1)-V(0,0) )*V(0,0) );
      
      SLocal(1,1)=  invvol * (  V(1,1)*V(1,1) + V(0,1)*V(0,1) );
      SLocal(1,2)=  invvol * ( -V(1,1)*V(1,0) - V(0,1)*V(0,0) );
      
      SLocal(2,2)=  invvol * ( V(1,0)*V(1,0)+ V(0,0)*V(0,0) );
      
      SLocal(1,0)=SLocal(0,1);
      SLocal(2,0)=SLocal(0,2);
      SLocal(2,1)=SLocal(1,2);
      
      // Quadrature for heat transfer coefficient
      double KLocal=(kappa(cells(icell,0))+kappa(cells(icell,1))+kappa(cells(icell,2)))/3.0;
      
      
      int i;
      for (int i=0;i<=ndim;i++)
        for (int j=0;j<=ndim;j++)
        {
          S(cells(icell,i),cells(icell,j))+=theta*KLocal*SLocal(i,j)+ tauinv*vol*MLocal0(i,j);
          Rhs(cells(icell,j))+=((1.0-theta)*KLocal*SLocal(i,j)+ tauinv*vol*MLocal0(i,j))* 
            OldSol(cells(icell,j))+vol*MLocal0(i,j)*source(cells(icell,j));
        }    
    }

    // Assemble boundary conditions (Dirichlet penalty method)
    int nbfaces=grid.nbfaces();
    auto bfaces=grid.get_bfaces();
    auto bfaceregions=grid.get_bfaceregions();
    
    for (int ibface=0; ibface<nbfaces; ibface++)
    {
      // Obtain number of boundary condition
      int ireg=bfaceregions(ibface);
      
      // Check if it is "Dirichlet"
      if (bcfac(ireg)>=Dirichlet)
      {
        // Assemble penalty values
        int i;
        i=bfaces(ibface,0);
        S(i,i)+=bcfac(ireg);
        Rhs(i)+=Dirichlet*bcval(ireg);
        
        i=bfaces(ibface,1);
        S(i,i)+=bcfac(ireg);
        Rhs(i)+=Dirichlet*bcval(ireg);
        
      }
      else if (bcfac(ireg)>0.0)
      {

        double dx= points(bfaces(ibface,1),0)-points(bfaces(ibface,0),0);
        double dy= points(bfaces(ibface,1),1)-points(bfaces(ibface,0),1);
        double h=sqrt(dx*dx+dy*dy);
        // first oder quadrature, needs to be replaced by second order
        int i0=bfaces(ibface,0);
        int i1=bfaces(ibface,1);
        
        S(i0,i0)+=theta*h*bcfac(ireg)/3.0;
        S(i0,i1)+=theta*h*bcfac(ireg)/6.0;
        
        S(i1,i0)+=theta*h*bcfac(ireg)/6.0;
        S(i1,i1)+=theta*h*bcfac(ireg)/3.0;
        
        Rhs(i0)+=OldSol(i0)*(1.0-theta)*h*bcfac(ireg)/3.0;
        Rhs(i0)+=OldSol(i1)*(1.0-theta)*h*bcfac(ireg)/6.0;
        Rhs(i1)+=OldSol(i0)*(1.0-theta)*h*bcfac(ireg)/6.0;
        Rhs(i1)+=OldSol(i1)*(1.0-theta)*h*bcfac(ireg)/3.0;
        
        Rhs(i0)+=0.5*h*bcval(ireg);
        Rhs(i0)+=0.5*h*bcval(ireg);
      }
    }
    S.flush();
  }

  double l2norm(std::shared_ptr<numcxx::SimpleGrid> pgrid, 
                std::shared_ptr<numcxx::DArray1> pu)
  {
    auto &grid=*pgrid;
    auto &u=*pu;
    auto ndim=grid.spacedim();
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map
    int npoints=grid.npoints();
    int ncells=grid.ncells();
    
    // Local mass matrix
    auto pMLocal0=numcxx::DMatrix::create({{2,1,1},{1,2,1},{1,1,2}}); // from Stroud quadrature...
    auto &MLocal0=*pMLocal0;
    MLocal0*=1.0/12.0;
    
    double norm=0.0;
    for (int icell=0; icell<ncells; icell++)
    {
      double vol;
      compute_cell_volume(icell,points,cells,vol);
      
      for (int i=0;i<=ndim;i++)
        for (int j=0;j<=ndim;j++)
          norm+=u(cells(icell,i))*u(cells(icell,j))*vol*MLocal0(i,j);
    }
    return sqrt(norm);
  }

  double h1norm(std::shared_ptr<numcxx::SimpleGrid> pgrid, 
                std::shared_ptr<numcxx::DArray1> pu)
  {
    auto &grid=*pgrid;
    auto &u=*pu;
    auto ndim=grid.spacedim();
    auto points=grid.get_points(); // Array of global nodes
    auto cells=grid.get_cells();   // Local-global dof map
    int npoints=grid.npoints();
    int ncells=grid.ncells();


    auto pSLocal=numcxx::DArray2::create(ndim+1, ndim+1);
    auto &SLocal=*pSLocal;

    double norm=0.0;
    for (int icell=0; icell<ncells; icell++)
    {
      double vol;
      compute_local_stiffness_matrix(icell, points,cells, SLocal, vol);
      for (int i=0;i<=ndim;i++)
        for (int j=0;j<=ndim;j++)
          norm+=u(cells(icell,i))*u(cells(icell,j))*SLocal(i,j);
    }
    return sqrt(norm);
  }



////////////////////////////////////////////////////////

  





}
#endif

