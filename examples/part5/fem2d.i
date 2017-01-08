
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
}

