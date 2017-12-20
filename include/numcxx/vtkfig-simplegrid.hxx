/// 
/// \file vtkfig-simplegrid.hxx
///
/// Header for adapter beteween vtkfig dataset and simple grid
/// 
#include "vtkfigDataSet.h"
#include "numcxx/simplegrid.hxx"

namespace numcxx
{

  /// Create vtkfig dataset from simple grid
  #ifdef VTKFIG

  inline std::shared_ptr<vtkfig::DataSet> vtkfigDataSet(const SimpleGrid &Grid)
  {
    auto griddata=vtkfig::DataSet::New();
    griddata->SetSimplexGrid(2,Grid.get_points(), Grid.get_cells());
    griddata->SetSimplexGridBoundaryCells(Grid.get_bfaces());
    griddata->SetCellRegions(Grid.get_cellregions());
    griddata->SetBoundaryCellRegions(Grid.get_bfaceregions());
    return griddata;
  }

  inline std::shared_ptr<vtkfig::DataSet> vtkfigDataSet(const std::shared_ptr<SimpleGrid> pGrid)
  {
    return vtkfigDataSet(*pGrid);
  }
  #endif
  
  
}
