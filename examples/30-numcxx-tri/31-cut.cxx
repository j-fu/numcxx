///
/// \example 31-cut.cxx
///
///  Create a triangular grid with multiple regions and visualize it using vtkfig
/// 
#include <cstdio>
#include <iostream>
#include <ctime>
#include <cmath>
#include "numcxx/numcxx.hxx"
#include "numcxx/simplegrid.hxx"
#include "numcxx/vtkfig-simplegrid.hxx"

#include "vtkfig/vtkfigFrame.h"
#include "vtkfig/vtkfigDataSet.h"
#include "vtkfig/vtkfigGridView.h"



int main(void)
{
  auto pGeometry=numcxx::Geometry::create();
  pGeometry->set_points({
      {-2,0},
      {0,0},
      {2,0},
      {2,2},
      {0.5,2},
      {0,1},
      {-0.5,2},
      {-2,2}
    });
  pGeometry->set_bfaces({
      {0,1},
      {1,2},
      {2,3},
      {3,4},
      {4,5},
      {5,6},
      {6,7},
      {7,0},
      {5,1},
        });
  
  
  pGeometry->set_bfaceregions({1,1,2,3,4,4,3,2,5});
  
  pGeometry->set_regionpoints({
      {-0.5,1},
      {0.5,1}
    });
  pGeometry->set_regionnumbers({1,2});
  pGeometry->set_regionvolumes({0.1,0.01});
  
  
  
  auto pGrid=numcxx::SimpleGrid::create(pGeometry,"zpaAqV");
  
  
#ifdef VTKFIG
  auto griddata=numcxx::vtkfigDataSet(pGrid);
  auto frame=vtkfig::Frame::New();
  auto gridview=vtkfig::GridView::New();
  gridview->SetData(griddata);
  frame->AddFigure(gridview);
  frame->Interact();
#endif
}


