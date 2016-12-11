%module mymodule

%{
#define SWIG_FILE_WITH_INIT

#include "numcxx/numcxx.hxx"
#include "numcxx/simplegrid.hxx"
#include "mymodule.h"

%}

template<class A> struct std::shared_ptr
{
    A* operator->() const;
};

namespace mymodule
{
   inline std::shared_ptr<numcxx::DArray1>  myfunction(std::shared_ptr<numcxx::SimpleGrid> g);
}

