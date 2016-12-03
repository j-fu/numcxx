#include <cstring>
#include <iostream>

namespace triangle
{
    extern "C" { 
#define VOID void
#include "../../triangle/triangle.h"
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
    
    inline SimpleGrid::SimpleGrid(const Geometry & geometry, const char * flags)
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
            int i;
            nregs=nholes=0;
            if (!geometry.regionnumbers)
            {
                throw std::runtime_error("numcxx: triangulate: Missing region numbers"); 
            }
            if (!geometry.regionvolumes)
            {
                throw std::runtime_error("numcxx: triangulate: Missing region volumes"); 
            }
            for (i=0;i<geometry.regionnumbers->shape(0);i++)
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
            
            for (i=0;i<geometry.regionnumbers->shape(0);i++)
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
            int *triattr,i;
            
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
}
