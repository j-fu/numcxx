/// 
/// \file vtkfig-simplegrid.hxx
///
/// Header for adapter beteween vtkfig dataset and simple grid
/// 
#include "vtkfig/vtkfigDataSet.h"
#include "numcxx/simplegrid.hxx"

namespace numcxx
{

  /// Create vtkfig dataset from simple grid
  #ifdef VTKFIG
  inline std::shared_ptr<vtkfig::DataSet> vtkfigDataSet(std::shared_ptr<SimpleGrid> pGrid)
  {
    auto griddata=vtkfig::DataSet::New();
    griddata->SetSimplexGrid(2,pGrid->get_points(), pGrid->get_cells());
    griddata->SetSimplexGridBoundaryCells(pGrid->get_bfaces());
    griddata->SetCellRegions(pGrid->get_cellregions());
    griddata->SetBoundaryCellRegions(pGrid->get_bfaceregions());
    return griddata;
  }
  #endif
  
  
}
