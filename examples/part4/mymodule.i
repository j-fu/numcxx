
/* Set the module name to be seen in python*/
%module mymodule

/* Add headers to be included from within the wrapped code*/
%{

#include "numcxx/numcxx.hxx"
#include "numcxx/simplegrid.hxx"
#include "mymodule.hxx"

%}


/* These are the actual declarations to be wrapped */
namespace mymodule
{
   inline std::shared_ptr<numcxx::DArray1>  myfunction(std::shared_ptr<numcxx::SimpleGrid> g);
}

