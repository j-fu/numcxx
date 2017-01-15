
/* Set the module name to be seen in python*/
%module fem2d

/* Add headers to be included from within the wrapped code*/
%{

#include "numcxx/numcxx.hxx"
#include "numcxx/simplegrid.hxx"
#include "fem2d.hxx"

%}


/* These are the actual declarations to be wrapped */
namespace fem2d
{
    inline std::shared_ptr<numcxx::DSparseMatrix> assemble_heat_matrix(std::shared_ptr<numcxx::SimpleGrid> g, std::shared_ptr<numcxx::IArray1> diri);
    inline std::shared_ptr<numcxx::DArray1> assemble_heat_rhs(std::shared_ptr<numcxx::SimpleGrid> g, std::shared_ptr<numcxx::IArray1> diri, std::shared_ptr<numcxx::DArray1> bcval);

  inline std::shared_ptr<numcxx::DSparseMatrix> 
  assemble_general_heat_matrix(
    std::shared_ptr<numcxx::SimpleGrid> grid,  // Discretization grid
    std::shared_ptr<numcxx::DArray1> kappa, // heat conduction coefficient (per node)
    std::shared_ptr<numcxx::DArray1> alpha // boundary heat transfer coefficient (per boundary region, value >=DirichletPenalty marks Dirichlet)
    );

  const double DirichletPenalty=1.0e30;

  // Assemble zero right hand side for mixed Dirichlet/Homogeneus Neumann problem
  std::shared_ptr<numcxx::DArray1>
  assemble_general_heat_rhs(
    std::shared_ptr<numcxx::SimpleGrid> grid,  // Discretization grid
    std::shared_ptr<numcxx::DArray1> f,        // heat source (per node)
    std::shared_ptr<numcxx::DArray1> g,        // boundary ambient temperature
    std::shared_ptr<numcxx::DArray1> alpha     // boundary heat transfer coefficient (large value marks Dirichlet)
    );

  void assemble_transient_heat_matrix_and_rhs(
    std::shared_ptr<numcxx::SimpleGrid> grid,// Discretization grid
    std::shared_ptr<numcxx::DSparseMatrix> S,
    std::shared_ptr<numcxx::DArray1> Rhs,
    std::shared_ptr<numcxx::DArray1> OldSol,
    std::shared_ptr<numcxx::DArray1> f,    // heat source (per node)
    std::shared_ptr<numcxx::DArray1> g,    // boundary ambient temperature
    std::shared_ptr<numcxx::DArray1> kappa, // heat conduction coefficient (per node)
    std::shared_ptr<numcxx::DArray1> alpha, // boundary heat transfer coefficient (per boundary region, value >=DirichletPenalty marks Dirichlet)
    double tau,
    double theta
    );



}

