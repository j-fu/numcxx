#ifndef NUMCXX_SIMPLEGRID_HXX
#define NUMCXX_SIMPLEGRID_HXX


#include "tarray1.hxx"
#include "tarray2.hxx"
#include "geometry.hxx"

namespace numcxx
{
  class SimpleGrid
  {
  public:
    SimpleGrid(){};
    SimpleGrid(const SimpleGrid& g)=delete;

    SimpleGrid(const Geometry &geometry, const char *triangle_flags);
    static std::shared_ptr<SimpleGrid> create(const Geometry &geometry, const std::string triangle_flags) { return std::make_shared<SimpleGrid>(geometry,triangle_flags.c_str());}
    static std::shared_ptr<SimpleGrid> create(std::shared_ptr<Geometry> geometry, const char * triangle_flags) { return std::make_shared<SimpleGrid>(*geometry,triangle_flags);}

    std::shared_ptr<TArray2<double>> points=nullptr;
    std::shared_ptr<TArray2<int>> cells=nullptr;
    std::shared_ptr<TArray1<int>> cellregions=nullptr;
    std::shared_ptr<TArray2<int>> bfaces=nullptr;
    std::shared_ptr<TArray1<int>> bfaceregions=nullptr;
        
    const TArray2<double>& get_points() { return *points;};
    const TArray2<int>& get_cells() { return *cells;};
    const TArray1<int>& get_cellregions() {return  *cellregions;};
    const TArray2<int>& get_bfaces() { return *bfaces;};
    const TArray1<int>& get_bfaceregions() { return *bfaceregions;};

    int spacedim() { return points->shape(1);}
    int griddim() { return bfaces->shape(1);}
    int ncells() {return cells->shape(0);}
    int npoints() {return points->shape(0);}
    int nbfaces() {return bfaces->shape(0);}

  };
}


#endif
