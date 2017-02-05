
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
      const double Dirichlet=1.0e30;
      
  inline void
  assemble_heat_problem_with_source(
    std::shared_ptr<numcxx::SimpleGrid> &pGrid,
    std::shared_ptr<numcxx::DArray1> pBCfac,
    std::shared_ptr<numcxx::DArray1> pBCval,
    std::shared_ptr<numcxx::DArray1> pSource,
    std::shared_ptr<numcxx::DArray1> pKappa,
    std::shared_ptr<numcxx::DSparseMatrix>&pS,
    std::shared_ptr<numcxx::DArray1> &pRhs
    );

  void assemble_and_apply_nonlinear_heat(
    std::shared_ptr<numcxx::SimpleGrid >pGrid,// Discretization grid
    std::shared_ptr<numcxx::DArray1> pBCfac,
    std::shared_ptr<numcxx::DArray1> pBCval,
    std::shared_ptr<numcxx::DArray1> pSol,
    std::shared_ptr<numcxx::DSparseMatrix> pS,
    std::shared_ptr<numcxx::DArray1> pRhs
    );

  void initialize_bc(
    std::shared_ptr<numcxx::SimpleGrid >grid,// Discretization grid
    std::shared_ptr<numcxx::DArray1> bcval,
    std::shared_ptr<numcxx::DArray1> Sol
    );

}

