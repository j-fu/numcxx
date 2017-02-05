
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

  const double Dirichlet=1.0e30;

  inline void
  assemble_simple_heat_problem(
    std::shared_ptr<numcxx::SimpleGrid> &pGrid,
    std::shared_ptr<numcxx::DArray1> pBCfac,
    std::shared_ptr<numcxx::DArray1> pBCval,
    std::shared_ptr<numcxx::DSparseMatrix>&pS,
    std::shared_ptr<numcxx::DArray1> &pRhs
    );

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

  void assemble_transient_heat_matrix_and_rhs(
    std::shared_ptr<numcxx::SimpleGrid> &pGrid,
    std::shared_ptr<numcxx::DArray1> pBCfac,
    std::shared_ptr<numcxx::DArray1> pBCval,
    std::shared_ptr<numcxx::DArray1> pSource,
    std::shared_ptr<numcxx::DArray1> pKappa,
    double tau, // time step
    double theta, //  choice of method
      bool lump,
    std::shared_ptr<numcxx::DArray1> pOldSol,
    std::shared_ptr<numcxx::DSparseMatrix>&pS,
    std::shared_ptr<numcxx::DArray1> &pRhs
    );


  double l2norm(std::shared_ptr<numcxx::SimpleGrid> pgrid, 
                std::shared_ptr<numcxx::DArray1> pu);

  double h1norm(std::shared_ptr<numcxx::SimpleGrid> pgrid, 
                std::shared_ptr<numcxx::DArray1> pu);


}

