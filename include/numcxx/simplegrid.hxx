///
/// \file simplegrid.hxx
///
/// Header for simple grid data class
/// 
#ifndef NUMCXX_SIMPLEGRID_HXX
#define NUMCXX_SIMPLEGRID_HXX


#include "tarray1.hxx"
#include "tarray2.hxx"
#include "geometry.hxx"

namespace numcxx
{
  ///
  /// Class containing data for simple grid data structure
  /// 
  class SimpleGrid
  {
  public:

    /// Trivial construtor
    SimpleGrid(){};
    SimpleGrid(const SimpleGrid& g)=delete;

    /// 
    /// Construct simple grid from geometry via library call to Triangle mesh generator
    /// 
    /// The \p triangle_flags \p must contain the z flag
    ///
    SimpleGrid(const Geometry &geometry, const char *triangle_flags);

    /// Static constructor from reference to geometry
    static std::shared_ptr<SimpleGrid> create(const Geometry &geometry, const std::string triangle_flags) { return std::make_shared<SimpleGrid>(geometry,triangle_flags.c_str());}

    /// Static constructor from shared pointer to geometry
    static std::shared_ptr<SimpleGrid> create(std::shared_ptr<Geometry> geometry, const char * triangle_flags) { return std::make_shared<SimpleGrid>(*geometry,triangle_flags);}


    /// Get array of point coordinates
    const TArray2<double>& get_points() const { return *points;}; 
    
    /// Get array of point indices describing cells
    const TArray2<int>& get_cells() const { return *cells;};

    /// Get array region markers
    const TArray1<int>& get_cellregions() const {return  *cellregions;};

    /// Get array of point indices describing boundary faces
    const TArray2<int>& get_bfaces() const { return *bfaces;};

    /// Get array of boundary markers
    const TArray1<int>& get_bfaceregions() const { return *bfaceregions;};

    /// Return dimension of space
    const int spacedim()  const { return points->shape(1);}

    /// Return dimension of grid
    const int griddim() const { return bfaces->shape(1);}

    /// Return number of cells
    const int ncells() const {return cells->shape(0);}

    /// Return number of points
    const int npoints() const {return points->shape(0);}

    /// Return number of boundary faces
    const int nbfaces() const {return bfaces->shape(0);}
    
  private:
    std::shared_ptr<TArray2<double>> points=nullptr;
    std::shared_ptr<TArray2<int>> cells=nullptr;
    std::shared_ptr<TArray1<int>> cellregions=nullptr;
    std::shared_ptr<TArray2<int>> bfaces=nullptr;
    std::shared_ptr<TArray1<int>> bfaceregions=nullptr;


  };
}


#endif
