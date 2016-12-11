#ifndef MYMODULE_H
#define MYMODULE_H

#include <numcxx/simplegrid.hxx>
#include <cmath>

namespace mymodule
{
    // This could be a function which solves a finite element problem
    // on the given grid
    // The c++11 convention says: if we do not want to store the grid in an object
    // created within this function, pass the grid as a reference
    inline std::shared_ptr<numcxx::TArray1<double>>  myfunction(numcxx::SimpleGrid &g)
    {
        
        auto points=g.get_points();
        int npoints=g.npoints();
        auto pF=numcxx::TArray1<double>::create(npoints);
        auto & F=*pF;
        for (int i=0;i<npoints; i++)
        {
            F(i)=std::sin(3.0*points(i,0))*std::cos(3.0*points(i,1));
        }
        return pF;
    }
    
    // Interfaceing with python requires to work with smart pointers... 
    inline std::shared_ptr<numcxx::TArray1<double>>  myfunction(std::shared_ptr<numcxx::SimpleGrid> pg) {return myfunction(*pg);}
    
}
#endif
