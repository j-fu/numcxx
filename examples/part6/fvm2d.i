
/* Set the module name to be seen in python*/
%module fvm2d

/* Add headers to be included from within the wrapped code*/
%{

#include "numcxx/numcxx.hxx"
#include "numcxx/simplegrid.hxx"
#include "fvm2d.hxx"

%}


/* These are the actual declarations to be wrapped */
namespace fvm2d
{
      

  inline std::shared_ptr<numcxx::DSparseMatrix> 
  assemble_general_heat_matrix(
    std::shared_ptr<numcxx::SimpleGrid> grid,  // Discretization grid
    std::shared_ptr<numcxx::DArray1> kappa, // heat conduction coefficient (per node)
    std::shared_ptr<numcxx::DArray1> alpha // boundary heat transfer coefficient (per boundary region, value >=DirichletPenalty marks Dirichlet)
    );

  const double DirichletPenalty=1.0e30;

  // Assemble zero right hand side for mixed Dirichlet/Homogeneus Neumann problem
  std::shared_ptr<numcxx::DArray1>
  assemble_heat_rhs(
    std::shared_ptr<numcxx::SimpleGrid> grid,  // Discretization grid
    std::shared_ptr<numcxx::DArray1> alpha,     // boundary heat transfer coefficient (large value marks Dirichlet)
    std::shared_ptr<numcxx::DArray1> g     // boundary temperature
    );


  void assemble_and_apply_nonlinear_heat(
    std::shared_ptr<numcxx::SimpleGrid >grid,// Discretization grid
    std::shared_ptr<numcxx::DSparseMatrix> S,
    std::shared_ptr<numcxx::DArray1> Sol,
    std::shared_ptr<numcxx::DArray1> Rhs,
    std::shared_ptr<numcxx::DArray1> alpha,
    std::shared_ptr<numcxx::DArray1> g
    );

  void initialize_bc(
    std::shared_ptr<numcxx::SimpleGrid >grid,// Discretization grid
    std::shared_ptr<numcxx::DArray1> Sol,
    std::shared_ptr<numcxx::DArray1> g
    );

}

