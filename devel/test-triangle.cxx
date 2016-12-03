#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include "numcxx/numcxx.hxx"



int main(void)
{
    auto pGeometry=numcxx::Geometry::create();
    pGeometry->set_points({
            {0,0},
            {1,0},
            {1,1}});
    pGeometry->set_bfaces({
            {0,1},
            {1,2},
            {2,0}});
    pGeometry->set_bfaceregions({0,1,2});

    pGeometry->set_regionpoints({{0.5,0.5}});
    pGeometry->set_regionnumbers({1});
    pGeometry->set_regionvolumes({0.01});

    
        
    auto pGrid=numcxx::SimpleGrid::create(pGeometry,"zpaAqV");
    std::cout<< *pGrid->bfaceregions;
}


