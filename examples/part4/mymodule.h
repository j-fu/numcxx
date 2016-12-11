#ifndef MYMODULE_H
#define MYMODULE_H

#include <numcxx/simplegrid.hxx>
#include <cmath>

namespace mymodule
{
inline std::shared_ptr<numcxx::TArray1<double>>  myfunction(numcxx::SimpleGrid &g)
{
    
    std::cout<< "myfunction\n";
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

inline std::shared_ptr<numcxx::TArray1<double>>  myfunction(std::shared_ptr<numcxx::SimpleGrid> pg) {return myfunction(*pg);}

}
#endif
