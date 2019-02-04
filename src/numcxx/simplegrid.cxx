#include <cstring>
#include <iostream>
#include <numcxx/numcxx.hxx>
#include <numcxx/simplegrid.hxx>

namespace triangle
{
  extern "C" { 
#define VOID void
#define REAL double 
#define ANSI_DECLARATORS
#include "../triangle/triangle.h"
#undef REAL
#undef VOID
#undef ANSI_DECLARATORS
  }
}

namespace numcxx
{
  inline void print_triangulateio(triangle::triangulateio *in)
  {
    std::printf("pointlist=%p\n",in->pointlist);
    std::printf("pointattributelist=%p\n",in->pointattributelist);
    std::printf("pointmarkerlist=%p\n",in->pointmarkerlist);
    std::printf("numberofpoints=%d\n",in->numberofpoints);
    std::printf("numberofpointattributes=%d\n",in->numberofpointattributes);
    std::printf("trianglelist=%p\n",in->trianglelist);
    std::printf("triangleattributelist=%p\n",in->triangleattributelist);
    std::printf("trianglearealist=%p\n",in->trianglearealist);
    std::printf("neighborlist=%p\n",in->neighborlist);
    std::printf("numberoftriangles=%d\n",in->numberoftriangles);
    std::printf("numberofcorners=%d\n",in->numberofcorners);
    std::printf("numberoftriangleattributes=%d\n",in->numberoftriangleattributes);
    std::printf("segmentlist=%p\n",in->segmentlist);
    std::printf("segmentmarkerlist=%p\n",in->segmentmarkerlist);
    std::printf("numberofsegments=%d\n",in->numberofsegments);
    std::printf("holelist=%p\n",in->holelist);
    std::printf("numberofholes=%d\n",in->numberofholes);
    std::printf("regionlist=%p\n",in->regionlist);
    std::printf("numberofregions=%d\n",in->numberofregions);
    std::printf("edgelist=%p\n",in->edgelist);
    std::printf("edgemarkerlist=%p\n",in->edgemarkerlist);
    std::printf("normlist=%p\n",in->normlist);
    std::printf("numberofedges=%d\n",in->numberofedges);
  }

  inline double dist(const numcxx::DArray2&points, int p1, int p2)
  {
    double dx=points(p2,0)-points(p1,0);
    double dy=points(p2,1)-points(p1,1);
    return sqrt(dx*dx+dy*dy);
  }

  void SimpleGrid::calc_hminmax(double& hmin, double& hmax) const
  {
    hmax=0.0;
    hmin=1.0e100;
    auto cells=*SimpleGrid::cells;
    auto points=*SimpleGrid::points;

    for (int icell=0; icell<ncells();icell++)
    {
      double h=dist(points, cells(icell,0), cells(icell,1));
      hmin=std::min(h,hmin);
      hmax=std::max(h,hmax);
   
      h=dist(points, cells(icell,1), cells(icell,2));
      hmin=std::min(h,hmin);
      hmax=std::max(h,hmax);

      h=dist(points, cells(icell,0), cells(icell,2));
      hmin=std::min(h,hmin);
      hmax=std::max(h,hmax);
    }

  }
    
  SimpleGrid::SimpleGrid(const Geometry & geometry, const char * flags)
  {
    if (strchr(flags,'z')==nullptr)
    {
      throw std::runtime_error("numcxx: triangulate: Missing 'z' flag"); 
    }
    
    struct triangle::triangulateio in,out;
    int dbg=0;
        
    out.pointlist=0;
    out.pointattributelist=0;
    out.pointmarkerlist=0;
    out.numberofpoints=0;
    out.numberofpointattributes=0;
    out.trianglelist=0;
    out.triangleattributelist=0;
    out.trianglearealist=0;
    out.neighborlist=0;
    out.numberoftriangles=0;
    out.numberofcorners=0;
    out.numberoftriangleattributes=0;
    out.segmentlist=0;
    out.segmentmarkerlist=0;
    out.numberofsegments=0;
    out.holelist=0;
    out.numberofholes=0;
    out.regionlist=0;
    out.numberofregions=0;
    out.edgelist=0;
    out.edgemarkerlist=0;
    out.normlist=0;
    out.numberofedges=0;
        
    in.pointlist=0;
    in.pointattributelist=0;
    in.pointmarkerlist=0;
    in.numberofpoints=0;
    in.numberofpointattributes=0;
    in.trianglelist=0;
    in.triangleattributelist=0;
    in.trianglearealist=0;
    in.neighborlist=0;
    in.numberoftriangles=0;
    in.numberofcorners=0;
    in.numberoftriangleattributes=0;
    in.segmentlist=0;
    in.segmentmarkerlist=0;
    in.numberofsegments=0;
    in.holelist=0;
    in.numberofholes=0;
    in.regionlist=0;
    in.numberofregions=0;
    in.edgelist=0;
    in.edgemarkerlist=0;
    in.normlist=0;
    in.numberofedges=0;
        
        
    // take  point list
    if (geometry.points)
    {
      in.pointlist=geometry.points->data();
      in.numberofpoints=geometry.points->shape(0);
      in.numberofpointattributes=0;
    }
    else
    {
      throw std::runtime_error("numcxx: triangulate: Missing point list"); 
    }
    // take bfaces if defined
    if (geometry.bfaces)
    {
      if (!geometry.bfaceregions)
        throw std::runtime_error("numcxx: triangulate: Missing bfaceregions"); 

      in.segmentlist=geometry.bfaces->data();
      in.segmentmarkerlist=geometry.bfaceregions->data();
      in.numberofsegments=geometry.bfaces->shape(0);
    }
    else
    {
      throw std::runtime_error("numcxx: triangulate: Missing bface list"); 
    }
        
    // create region and hole information for triangle
    // Holes are marked with region number <0
    if (geometry.regionpoints)
    {
      int nregs,nholes,ireg,ihole;
      nregs=nholes=0;
      if (!geometry.regionnumbers)
      {
        throw std::runtime_error("numcxx: triangulate: Missing region numbers"); 
      }
      if (!geometry.regionvolumes)
      {
        throw std::runtime_error("numcxx: triangulate: Missing region volumes"); 
      }
      for (long i=0;i<geometry.regionnumbers->shape(0);i++)
      {
        if ((*geometry.regionnumbers)(i)>=0)
          nregs++; 
        else 
          nholes++;
      }
      in.numberofholes=nholes;
      in.numberofregions=nregs;
      if(nholes)
        in.holelist=(double*)malloc(2*nholes*sizeof(double));
      in.regionlist=(double*)malloc(4*nregs*sizeof(double));
      ireg=ihole=0;
            
      for (long i=0;i<geometry.regionnumbers->shape(0);i++)
        if ((*geometry.regionnumbers)(i)>=0)
        {
          in.regionlist[ireg+0]=(*geometry.regionpoints)(i,0);
          in.regionlist[ireg+1]=(*geometry.regionpoints)(i,1);
          in.regionlist[ireg+2]=(*geometry.regionnumbers)(i);
          in.regionlist[ireg+3]=(*geometry.regionvolumes)(i);
          ireg+=4;
        }		 
        else
        {
          in.holelist[ihole+0]=(*geometry.regionpoints)(i,0); 
          in.holelist[ihole+1]=(*geometry.regionpoints)(i,1); 
          ihole+=2;
        }
            
    }    
        
    
    if (dbg)
    {
      printf("in:\n");
      print_triangulateio(&in);
    }
        
        
    //  Do the damn thing!
    triangle::triangulate((char*)flags,&in,&out,NULL);
        
    if (dbg)
    {
      printf("out (flags=%s):\n",flags);
      print_triangulateio(&out);
    }
        
    // This triangle specific input information is not needed anymore.
    if(in.holelist) triangle::trifree(in.holelist);
    if(in.regionlist) triangle::trifree(in.regionlist);
        
    
    // Take point and cell list as is. 
    points=std::make_shared<numcxx::TArray2<double>>(out.numberofpoints,2,out.pointlist,[](double*p){triangle::trifree(p);});
    cells=std::make_shared<numcxx::TArray2<int>>(out.numberoftriangles,3,out.trianglelist,[](int*p){triangle::trifree(p);});
        

    // Our cellregions are integer, not double as in triangle. 
    // Furthermore, they are mandatory in our process, so we
    // set them to 0 if they are not available
    {
      int i;
            
      cellregions=std::make_shared<numcxx::TArray1<int>>(out.numberoftriangles);
      if (out.triangleattributelist)
        for (i=0;i<out.numberoftriangles;i++) 
          (*cellregions)(i)=(int)out.triangleattributelist[i];
      else
        for (i=0;i<out.numberoftriangles;i++) 
          (*cellregions)(i)=0;
    
      if (in.triangleattributelist) triangle::trifree(in.triangleattributelist);
      if (out.triangleattributelist) triangle::trifree(out.triangleattributelist);
    }
        
    if (out.segmentlist)
    {
      bfaces=std::make_shared<numcxx::TArray2<int>>(out.numberofsegments,2,out.segmentlist,[](int*p){triangle::trifree(p);});
      bfaceregions=std::make_shared<numcxx::TArray1<int>>(out.numberofsegments,out.segmentmarkerlist,[](int*p){triangle::trifree(p);});
    }
    else
    {
      throw std::runtime_error("numcxx: triangulate: Missing segment list"); 
    }

  }


  /// 
  /// Construct simple grid from array of x/y coordinats
  ///

  bool is_monotone(const DArray1 & x)
  {
    int nx=x.shape(0);
    for (int i=1;i<nx;i++)
    {
      if (x(i)<x(i-1)) return false;
    }
    return true;
  }

  inline bool leq(double x, double x1, double x2)
  {
    if (x>x1) return false;
    if (x>x2) return false;
    return true;
  }
  
  inline bool geq(double x, double x1, double x2)
  {
    if (x<x1) return false;
    if (x<x2) return false;
    return true;
  }

  inline void check_bface(int n1, int n2,
                          double x1,double xn,
                          double y1,double yn,
                          DArray2& coord, IArray2& BFN, IArray1 & BFR,  
                          int &ibf)
  {
    if (geq(x1,coord(n1,0),coord(n2,0))) { BFN(ibf,0)=n1; BFN(ibf,1)=n2; BFR(ibf)=4; ++ibf; return;}
    if (leq(xn,coord(n1,0),coord(n2,0))) { BFN(ibf,0)=n1; BFN(ibf,1)=n2; BFR(ibf)=2; ++ibf; return;}
    if (geq(y1,coord(n1,1),coord(n2,1))) { BFN(ibf,0)=n1; BFN(ibf,1)=n2; BFR(ibf)=1; ++ibf; return;}
    if (leq(yn,coord(n1,1),coord(n2,1))) { BFN(ibf,0)=n1; BFN(ibf,1)=n2; BFR(ibf)=3; ++ibf; return;}
  }                                                                                         
  
  
  
  SimpleGrid::SimpleGrid(const DArray1 & x, const  DArray1 & y)
  {
    int nx=x.shape(0);
    int ny=y.shape(0);
    assert(is_monotone(x));
    assert(is_monotone(y));
    
    
    double hmin=x(1)-x(0);
    for (int i=0;i<nx-1;i++) if ((x(i+1)-x(i)) <hmin) hmin=x(i+1)-x(i);
    for (int i=0;i<ny-1;i++) if ((y(i+1)-y(i)) <hmin) hmin=y(i+1)-y(i);
    
    assert(hmin>0.0);
    double eps=1.0e-5*hmin;
    
    int nnodes=nx*ny;
    int ncells=2*(nx-1)*(ny-1);
    int nbfaces=2*(nx-1)+2*(ny-1);
    
    points=DArray2::create(nnodes,2);
    cells=IArray2::create(ncells,3);
    cellregions=IArray1::create(ncells);
    bfaces=IArray2::create(nbfaces,2);
    bfaceregions=IArray1::create(nbfaces);
    
    
    auto & coord=*points;
    int icoord=0;
    for(int iy=0;iy<ny;iy++)
      for (int ix=0;ix<nx;ix++)
      {
        coord(icoord,0)=x(ix);
        coord(icoord,1)=y(iy);
        icoord++;
      }
    assert(icoord==nnodes);
    
    auto & cn=*cells;
    auto & cr=*cellregions;
    int icell=0;
    for(int iy=0;iy<ny-1;iy++)
      for (int ix=0;ix<nx-1;ix++)
      {
	int ip=ix+iy*nx;
	int p00 = ip;
	int p10 = ip+1;
        int p01 = ip  +nx;
        int p11 = ip+1+nx;
        
        cn(icell,0)=p00;
        cn(icell,1)=p10;
        cn(icell,2)=p11;
        cr(icell)=1;
        icell++;
        cn(icell,0)=p11;
        cn(icell,1)=p01;
        cn(icell,2)=p00;
        cr(icell)=1;
        icell++;
        
      }
    assert(icell==ncells);
    
    // lazy way to  create boundary grid
    double x1=x(0)+eps;
    double xn=x(nx-1)-eps;
    double y1=y(0)+eps;
    double yn=y(ny-1)-eps;
    
    int ibface=0;
    // lazy but easy...
    for (int icell=0;icell<ncells;icell++)
    {
      int n1=cn(icell,0);
      int n2=cn(icell,1);
      int n3=cn(icell,2);
      check_bface(n1,n2,x1,xn,y1,yn,coord,*bfaces,*bfaceregions,ibface);
      check_bface(n1,n3,x1,xn,y1,yn,coord,*bfaces,*bfaceregions,ibface);
      check_bface(n2,n3,x1,xn,y1,yn,coord,*bfaces,*bfaceregions,ibface);
    }
    assert(ibface==nbfaces);
  }
  
  
  
}




