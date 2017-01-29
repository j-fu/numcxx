
/* Set the module name to be seen in python*/
%module convdiff

/* Add headers to be included from within the wrapped code*/
%{

#include "numcxx/numcxx.hxx"
#include "numcxx/simplegrid.hxx"
#include "convdiff.hxx"

%}


/* These are the actual declarations to be wrapped */
namespace convdiff
{
      std::shared_ptr<numcxx::DArray1> solve_expfit(int n, double v, double D);
      std::shared_ptr<numcxx::DArray1> solve_central(int n, double v, double D);
      std::shared_ptr<numcxx::DArray1> solve_upwind(int n, double v, double D);

}
