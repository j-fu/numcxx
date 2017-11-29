///
/// \file geometry.hxx
///
/// Header for geometry description
/// 
#ifndef NUMCXX_GEOMETRY_H
#define NUMCXX_GEOMETRY_H

#include "tarray1.hxx"
#include "tarray2.hxx"

namespace numcxx
{
  ///
  /// Class collecting data for the description of piecewise linear geometries
  /// 
  class Geometry
  {
        
  public:

    /// Set member via intializer list
    void set_points       (const std::initializer_list<std::initializer_list<double>> &il) {points=TArray2<double>::create(il);};

    /// Set member via intializer list
    void set_bfaces       (const std::initializer_list<std::initializer_list<int>> &il) {bfaces=TArray2<int>::create(il);};

    /// Set member via intializer list
    void set_bfaceregions (const std::initializer_list<int> &il){bfaceregions=TArray1<int>::create(il);};

    /// Set member via intializer list
    void set_regionpoints (const std::initializer_list<std::initializer_list<double>> &il) {regionpoints=TArray2<double>::create(il);};

    /// Set member via intializer list
    void set_regionnumbers(const std::initializer_list<int> &il){regionnumbers=TArray1<int>::create(il);};

    /// Set member via intializer list
    void set_regionvolumes(const std::initializer_list<double> &il){regionvolumes=TArray1<double>::create(il);};


    ///
    /// Points: npt x dim array of double containing point coordinates
    ///
    std::shared_ptr<TArray2<double>> points=nullptr;

    ///
    /// nbfaces x dim array of  of integers describing boundary segments
    ///
    std::shared_ptr<TArray2<int>>    bfaces=nullptr;


    ///
    /// nbfaces array of integers describing  boundary segment markers
    ///
    std::shared_ptr<TArray1<int>>    bfaceregions=nullptr;

    ///
    /// nreg x dim array of doubles containing  point coordinates of region points
    ///
    std::shared_ptr<TArray2<double>> regionpoints=nullptr;

    ///
    /// nreg array of integers containing region markers
    ///
    std::shared_ptr<TArray1<int>>    regionnumbers=nullptr;

    ///
    /// nreg array of integers  containing the maximum volumes/areas of triangles in a region
    ///
    std::shared_ptr<TArray1<double>> regionvolumes=nullptr;


    static std::shared_ptr<Geometry> create() { return std::make_shared<Geometry>();}
    static std::shared_ptr<Geometry> New() { return std::make_shared<Geometry>();}
  };
    
}


#endif
