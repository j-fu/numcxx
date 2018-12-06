///
/// \example 30-triangle.cxx
///
///  Create a triangular grid and visualize it using vtkfig
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
            {0,0},
            {1,0},
            {1,1}});
    pGeometry->set_bfaces({
            {0,1},
            {1,2},
            {2,0}});
    pGeometry->set_bfaceregions({1,2,3});

    pGeometry->set_regionpoints({{0.5,0.5}});
    pGeometry->set_regionnumbers({1});
    pGeometry->set_regionvolumes({0.01});

    
        
    auto pGrid=numcxx::SimpleGrid::create(pGeometry,"zpaAqV");


    #ifdef VTKFIG
    auto griddata=numcxx::vtkfigDataSet(pGrid);
    auto frame=vtkfig::Frame::New();
    auto gridview=vtkfig::GridView::New();
    gridview->SetData(griddata);
    frame->AddFigure(gridview);
    frame->Show();
    frame->Interact();
    #endif
}


