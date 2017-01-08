/* -*- C -*-  (not really, but good for syntax highlighting) */
%module numcxxwrap

/*
Section 1: Helper functions to convert data between C++/numpy
*/
%{
#define SWIG_FILE_WITH_INIT
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numcxx/numcxx.hxx"
#include "numcxx/simplegrid.hxx"
    namespace numcxx
    {
        // Proxy class for numpy object to be stored 
        class NumpyProxy
        {
            PyObject *o;
        public:
            NumpyProxy(PyObject *o): o(o) { Py_INCREF(o);};
            ~NumpyProxy() {Py_DECREF(o);};
        };
        
        // Convert numpy object with data to numcxx
        std::shared_ptr<DArray1 > NumpyAsDArray1(PyObject *o,double* v, int n1)
        {
            auto proxy=std::make_shared<NumpyProxy>(o);
            return std::make_shared<DArray1 >(n1,v,proxy);
        }
        
        // Convert numcxx data to numpy
        void DArray1AsNumpy(std::shared_ptr<DArray1 > a,PyObject *o,double** ARGOUTVIEW_ARRAY1, int* DIM1)
        {
            *DIM1=a->size();
            *ARGOUTVIEW_ARRAY1=a->data();
        }

        std::shared_ptr<DArray2 > NumpyAsDArray2(PyObject *o,double* v, int n1, int n2)
        {
            auto proxy=std::make_shared<NumpyProxy>(o);
            return std::make_shared<DArray2>(n1, n2, v,proxy);
        }

        void DArray2AsNumpy(std::shared_ptr<DArray2 > a,PyObject *o,double** ARGOUTVIEW_ARRAY2, int* DIM1,int* DIM2)
        {
            *DIM1=a->shape(0);
            *DIM2=a->shape(1);
            *ARGOUTVIEW_ARRAY2=a->data();
        }


        std::shared_ptr<IArray1 > NumpyAsIArray1(PyObject *o,int* v, int n1)
        {
            auto proxy=std::make_shared<NumpyProxy>(o);
            return std::make_shared<IArray1 >(n1,v,proxy);
        }
        
        void IArray1AsNumpy(std::shared_ptr<IArray1 > a,PyObject *o,int** ARGOUTVIEW_ARRAY1, int* DIM1)
        {
            *DIM1=a->size();
            *ARGOUTVIEW_ARRAY1=a->data();
        }

        std::shared_ptr<IArray2 > NumpyAsIArray2(PyObject *o,int* v, int n1, int n2)
        {
            auto proxy=std::make_shared<NumpyProxy>(o);
            return std::make_shared<IArray2>(n1, n2, v,proxy);
        }

        void IArray2AsNumpy(std::shared_ptr<IArray2 > a,PyObject *o,int** ARGOUTVIEW_ARRAY2, int* DIM1,int* DIM2)
        {
            *DIM1=a->shape(0);
            *DIM2=a->shape(1);
            *ARGOUTVIEW_ARRAY2=a->data();
        }


        std::shared_ptr<DMatrix> NumpyAsDMatrix(PyObject *o,double* v, int n1, int n2)
        {
            auto proxy=std::make_shared<NumpyProxy>(o);
            return std::make_shared<DMatrix>(n1,n2, v,proxy);
        }

        void DMatrixAsNumpy(std::shared_ptr<DMatrix > a,PyObject *o,double** ARGOUTVIEW_ARRAY2, int* DIM1,int* DIM2)
        {
            *DIM1=a->shape(0);
            *DIM2=a->shape(1);
            *ARGOUTVIEW_ARRAY2=a->data();
        }
    }
    
%}


/*
  Section 2: include standard typemap files
 */
%include <std_vector.i>
%include <numpy.i>
%clear(double** ARGOUTVIEW_ARRAY1, int* DIM1);
%clear(int** ARGOUTVIEW_ARRAY1, int* DIM1);





/* Section 3: Typemaps for numcxx <-> numpy array communication.
   
   This typemap  file needs to be included after the original,
   unmodified numpy.i.

   The general idea is to catch  the python objects in order to 
   use Py_INCREF/Py_DECREF in the proper places.
   
*/

%define %numcxx_numpy_typemaps(DATA_TYPE, DATA_TYPECODE, DIM_TYPE)

 /***************************************************************/
 /***************************************************************/
 /***************************************************************/
 /*  
     Numpy array as numcxx array.
 */

 /***************************************************************
   Typemap suite for (PyObject *O, DATA_TYPE* INPLACE_ARRAY1, int DIM1)
   for numpy arrays as numcxx arrays
 */


%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY, fragment="NumPy_Macros")
(PyObject* O, DATA_TYPE* INPLACE_ARRAY1, int DIM1)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input), DATA_TYPECODE);
}

%typemap(in,fragment="NumPy_Fragments")
  (PyObject* O, DATA_TYPE* INPLACE_ARRAY1, int DIM1)
  (PyArrayObject* array=NULL, int i=1)
{
  $1=$input;
  array = obj_to_array_no_conversion($input, DATA_TYPECODE);
  if (!array || 
      !require_dimensions(array,1) || 
      !require_contiguous(array) || 
      !require_native(array)) SWIG_fail;
  $2 = (DATA_TYPE*) array_data(array);
  $3 = array_size(array,0);
}


/***************************************************************
   Typemap suite for (PyObject *O, DATA_TYPE* INPLACE_ARRAY1, int DIM1, int DIM2)
   for numpy arrays as numcxx arrays
*/
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY, fragment="NumPy_Macros")
(PyObject* O, DATA_TYPE* INPLACE_ARRAY2, int DIM1, int DIM2)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input), DATA_TYPECODE);
}

%typemap(in,fragment="NumPy_Fragments")
  (PyObject* O, DATA_TYPE* INPLACE_ARRAY2, int DIM1, int DIM2)
  (PyArrayObject* array=NULL, int i1=1, int i2=1)
{
  $1=$input;
  array = obj_to_array_no_conversion($input, DATA_TYPECODE);
  if (!array || 
      !require_dimensions(array,2) ||
      // pycall.jl will have problems here but rightly so...
      !require_contiguous(array) ||  
      !require_native(array)) SWIG_fail;
  $2 = (DATA_TYPE*) array_data(array);
  $3 = array_size(array,0);
  $4 = array_size(array,1);
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
/*  
    numcxx array as  numpy: 
    
    Concerning PyArray_New  instead of PyArray_SimpleNewFromData, see the source on
    https://github.com/numpy/numpy/blob/master/numpy/core/src/multiarray/ctors.c:

    The obj (last) parameter is increfed and stored, so we pass the python proxy
    of the numcxx array.
*/

/***************************************************************
   Typemap suite for (PyObject *O, DATA_TYPE** ARGOUTVIEW_ARRAY1, int* DIM1)
   for numcxx array as  numpy.
   
   We use obj0 for the python obj passed after observation of the generated
   code. This should be fixed and declared properly in the parameter list.
*/
%typemap(in,numinputs=0)
  (PyObject *O,DATA_TYPE** ARGOUTVIEW_ARRAY1, int* DIM1)
  ( PyObject *o_tmp=NULL,DATA_TYPE*  data_temp = NULL,int  dim_temp )
{
  $1 = o_tmp;
  $2 = &data_temp;
  $3 = &dim_temp;
}
%typemap(argout, fragment="NumPy_Backward_Compatibility")
  (PyObject *O,DATA_TYPE** ARGOUTVIEW_ARRAY1, int* DIM1)
{
  npy_intp dims[1] = { *$3 };
  PyObject* obj = PyArray_SimpleNewFromData(1, dims, DATA_TYPECODE,(void*)(*$2));
  PyArrayObject* array = (PyArrayObject*) obj;
  if (!array) SWIG_fail;
  Py_INCREF(obj0);
  PyArray_SetBaseObject(array,obj0);
  $result = SWIG_Python_AppendOutput($result,obj);
}


/***************************************************************
   Typemap suite for (PyObject *O, DATA_TYPE** ARGOUTVIEW_ARRAY2, int* DIM1, int *DIM2)
   for numcxx array as  numpy.
   
   We use obj0 for the python obj passed after observation of the generated
   code. This should be fixed and declared properly in the parameter list.
*/
%typemap(in,numinputs=0)
 (PyObject *O,DATA_TYPE** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2)
  ( PyObject *o_tmp=NULL,DATA_TYPE*  data_temp = NULL,int  dim1_temp ,int  dim2_temp )
{
  $1 = o_tmp;
  $2 = &data_temp;
  $3 = &dim1_temp;
  $4 = &dim2_temp;
}
%typemap(argout, fragment="NumPy_Backward_Compatibility")
  (PyObject *O,DATA_TYPE** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2)
{
  npy_intp dims[2] = { *$3, *$4 };
  PyObject* obj = PyArray_SimpleNewFromData(2, dims, DATA_TYPECODE,(void*)(*$2));
  PyArrayObject* array = (PyArrayObject*) obj;
  if (!array) SWIG_fail;
  Py_INCREF(obj0);
  PyArray_SetBaseObject(array,obj0);
  $result = SWIG_Python_AppendOutput($result,obj);
}

%enddef    /* %numcxx_numpy_typemaps() macro */

%numcxx_numpy_typemaps(double, NPY_DOUBLE, int)
%numcxx_numpy_typemaps(int,NPY_INT, int)



%init %{
  import_array();
%}

/*
Section 4: decription of classes to be wrapped.
Essentially a subset of C++ header file.
 */

template<class A> struct std::shared_ptr
{
    A* operator->() const;
};
  
/* numcxx class wrappers */

namespace numcxx
{
    typedef unsigned int index;


    class DArray1
    {
    public:
        const index ndim();
        index size();
        int shape(int idim);
        
        static std::shared_ptr<DArray1> create(index n0);
        double item(index i0);
        void itemset(index i0, double x);
        
        std::shared_ptr<DArray1> copy();
        double __getitem__(index i) const;
        void __setitem__(index i,double);
        bool is_matrix();
    };
    
    class DArray2
    {
    public:
        const index ndim();
        index size();
        int shape(int idim);
        
        static std::shared_ptr<DArray2> create(index n0, index n1);
        double item(index i0,index i1);
        void itemset(index i0,index i1, double x);
        
        std::shared_ptr<DArray2> copy();
        std::shared_ptr<DArray1>   __getitem__(int i);
        bool is_matrix();
    };


    class IArray1
    {
    public:
        const index ndim();
        index size();
        int shape(int idim);
        
        static std::shared_ptr<IArray1> create(index n0);
        int  item(index i0);
        void itemset(index i0,int x);
        
        std::shared_ptr<IArray1> copy();
        int __getitem__(index i) const;
        void __setitem__(index i,int);
        bool is_matrix();
    };
    
    class IArray2
    {
    public:
        const index ndim();
        index size();
        int shape(int idim);
        
        static std::shared_ptr<IArray2> create(index n0, index n1);
        double item(index i0,index i1);
        void itemset(index i0,index i1, int x);
        
        std::shared_ptr<IArray2> copy();
        std::shared_ptr<IArray1>   __getitem__(int i);
        bool is_matrix();
    };


    class DSolverLapackLU
    {
    public:
        static std::shared_ptr<DSolverLapackLU> create(const std::shared_ptr<numcxx::DMatrix> pMatrix);
        void update();
        void solve( std::shared_ptr<DArray1> & Sol,  const std::shared_ptr<DArray1> & Rhs) const;
    };



    class DMatrix
    {
    public:
        const index ndim();
        index size();
        int shape(int idim);
        
        static std::shared_ptr<DMatrix> create(index n1, index n2);
        double item(index i0,index i1);
        void itemset(index i0,index i1, double x);
        
        std::shared_ptr<DMatrix> copy();
        std::shared_ptr<DArray1>   __getitem__(int i);
        bool is_matrix();
        std::shared_ptr<DMatrix> calculate_inverse();
    };
    


    class DSparseMatrix
    {
    public:
      static std::shared_ptr<DSparseMatrix> create(index n1, index n2);
      void  flush();
      std::shared_ptr< DMatrix > copy_as_dense();
      index shape(int idim);
    };

    class DSolverUMFPACK
    {
    public:
      static std::shared_ptr<DSolverUMFPACK> create(const std::shared_ptr<DSparseMatrix> pA);
      void update();
      void solve( std::shared_ptr<DArray1> Sol,  const std::shared_ptr<DArray1> Rhs);
    };

    



    class Geometry
    {
        
    public:
        
        std::shared_ptr<DArray2> points;
        std::shared_ptr<IArray2>    bfaces;
        std::shared_ptr<IArray1>    bfaceregions;
        std::shared_ptr<DArray2> regionpoints;
        std::shared_ptr<IArray1>    regionnumbers;
        std::shared_ptr<DArray1> regionvolumes;
        static std::shared_ptr<Geometry> create();
    };


    class SimpleGrid
    {
    public:
        static std::shared_ptr<SimpleGrid> create(std::shared_ptr<Geometry> geometry, const char *triangle_flags);
        std::shared_ptr<DArray2> points;
        std::shared_ptr<IArray2> cells;
        std::shared_ptr<IArray1> cellregions;
        std::shared_ptr<IArray2> bfaces;
        std::shared_ptr<IArray1> bfaceregions;

        int spacedim();
        int griddim();
        int ncells();  
        int npoints(); 
        int nbfaces();
        
    };

    
    std::shared_ptr<DArray1> NumpyAsDArray1(PyObject *O, double* INPLACE_ARRAY1, int DIM1);
    std::shared_ptr<DArray2> NumpyAsDArray2(PyObject *O, double* INPLACE_ARRAY2, int DIM1, int DIM2);
    std::shared_ptr<IArray1> NumpyAsIArray1(PyObject *O, int* INPLACE_ARRAY1, int DIM1);
    std::shared_ptr<IArray2> NumpyAsIArray2(PyObject *O, int* INPLACE_ARRAY2, int DIM1, int DIM2);
    std::shared_ptr<DMatrix> NumpyAsDMatrix(PyObject *O, double* INPLACE_ARRAY2, int DIM1, int DIM2);
    
    void DArray1AsNumpy(std::shared_ptr<DArray1> a,PyObject *O,double** ARGOUTVIEW_ARRAY1, int* DIM1);
    void DArray2AsNumpy(std::shared_ptr<DArray2> a,PyObject *O,double** ARGOUTVIEW_ARRAY2, int* DIM1,  int* DIM2);
    void IArray1AsNumpy(std::shared_ptr<IArray1> a,PyObject *O,int** ARGOUTVIEW_ARRAY1, int* DIM1);
    void IArray2AsNumpy(std::shared_ptr<IArray2> a,PyObject *O,int** ARGOUTVIEW_ARRAY2, int* DIM1,  int* DIM2);
    void DMatrixAsNumpy(std::shared_ptr<DMatrix> a,PyObject *O,double** ARGOUTVIEW_ARRAY2, int* DIM1,  int* DIM2);
    
}

/*
Section 5: template expansions
 */

%template(shared_ptrDArray1) std::shared_ptr<numcxx::DArray1 >;
%template(shared_ptrDArray2) std::shared_ptr<numcxx::DArray2 >;
%template(shared_ptrIArray1) std::shared_ptr<numcxx::IArray1 >;
%template(shared_ptrIArray2) std::shared_ptr<numcxx::IArray2 >;
%template(shared_ptrDMatrix) std::shared_ptr<numcxx::DMatrix >;
%template(shared_ptrDSparseMatrix) std::shared_ptr<numcxx::DSparseMatrix >;
%template(shared_ptrDSolverUMFPACK) std::shared_ptr<numcxx::DSolverUMFPACK >;
%template(shared_ptrDSolverLapackLU) std::shared_ptr<numcxx::DSolverLapackLU >;
%template(shared_ptrSimpleGrid) std::shared_ptr<numcxx::SimpleGrid >;
%template(shared_ptrGeometry) std::shared_ptr<numcxx::Geometry >;


/*
Section 6: Python code to extend functionality
 */
%pythoncode{
import numpy

"Shortcuts to numpy <-> numcxx array conversions"

def asdarray(proto):
    if type(proto)==numpy.ndarray:
        a=proto
    elif type(proto)==list:
        a=numpy.array(proto,dtype=numpy.float64)

    else:
        raise TypeError("asnumcxx: wrong type of input -- must be numpy.ndarray or list.")

    if a.dtype== numpy.float64 and a.ndim==1:
        return NumpyAsDArray1(a)          
    elif a.dtype== numpy.float64 and a.ndim==2:
        return NumpyAsDArray2(a)           

    raise TypeError("asnumcxx: dtype %s not supported"%a.dtype)

def asiarray(proto):
    if type(proto)==numpy.ndarray:
        a=proto
    elif type(proto)==list:
        a= numpy.array(proto,dtype=numpy.intc)
    else:
        raise TypeError("asnumcxx: wrong type of input -- must be numpy.ndarray or list.")

    if a.dtype==numpy.intc and a.ndim==1:
        return NumpyAsIArray1(a)          
    elif a.dtype==numpy.intc and a.ndim==2:
        return NumpyAsIArray2(a)           

    raise TypeError("asnumcxx: dtype %s not supported"%a.dtype)

def asdmatrix(proto):
    if type(proto)==numpy.ndarray:
        a=proto
    elif type(proto)==list:
        a=numpy.array(proto,dtype=numpy.float64)
    else:
        raise TypeError("asnumcxx: wrong type of input -- must be numpy.ndarray or list.")

    if a.dtype==numpy.float64 and a.ndim==2:
        return NumpyAsDMatrix(a)           

    raise TypeError("asnumcxx: error")

    
def asnumpy(a):
    if isinstance(a,shared_ptrDArray1):
        return DArray1AsNumpy(a)
    elif isinstance(a,shared_ptrDArray2):
        return DArray2AsNumpy(a)
    elif isinstance(a,shared_ptrIArray1):
        return IArray1AsNumpy(a)
    elif isinstance(a,shared_ptrIArray2):
        return IArray2AsNumpy(a)
    elif isinstance(a,shared_ptrDMatrix):
        return DMatrixAsNumpy(a)

    raise RuntimeError("asnumpy: error")

"Extension of geometry class"

def Geometry_set_points(self,proto):
    self.points=asdarray(proto)

def Geometry_set_bfaces(self,proto):
    self.bfaces=asiarray(proto)

def Geometry_set_bfaceregions(self,proto):
    self.bfaceregions=asiarray(proto)

def Geometry_set_regionpoints(self,proto):
    self.regionpoints=asdarray(proto)

def Geometry_set_regionnumbers(self,proto):
    self.regionnumbers=asiarray(proto)

def Geometry_set_regionvolumes(self,proto):
    self.regionvolumes=asdarray(proto)


def Geometry_get_points(self):
    return asnumpy(self.points)

def Geometry_get_bfaces(self):
    return asnumpy(self.bfaces)

def Geometry_get_bfaceregions(self):
    return asnumpy(self.bfaceregions)

def Geometry_get_regionpoints(self):
    return asnumpy(self.regionpoints)

def Geometry_get_regionnumbers(self):
    return asnumpy(self.regionnumbers)

def Geometry_get_regionvolumes(self):
    return asnumpy(self.regionvolumes)


shared_ptrGeometry.set_points=Geometry_set_points
shared_ptrGeometry.set_bfaces=Geometry_set_bfaces
shared_ptrGeometry.set_bfaceregions=Geometry_set_bfaceregions
shared_ptrGeometry.set_regionpoints=Geometry_set_regionpoints
shared_ptrGeometry.set_regionnumbers=Geometry_set_regionnumbers
shared_ptrGeometry.set_regionvolumes=Geometry_set_regionvolumes

shared_ptrGeometry.get_points=Geometry_get_points
shared_ptrGeometry.get_bfaces=Geometry_get_bfaces
shared_ptrGeometry.get_bfaceregions=Geometry_get_bfaceregions
shared_ptrGeometry.get_regionpoints=Geometry_get_regionpoints
shared_ptrGeometry.get_regionnumbers=Geometry_get_regionnumbers
shared_ptrGeometry.get_regionvolumes=Geometry_get_regionvolumes


"Extension of SimpleGrid class"

def SimpleGrid_get_points(self):
    return asnumpy(self.points)
def SimpleGrid_get_cells(self):
    return asnumpy(self.cells)
def SimpleGrid_get_cellregions(self):
    return asnumpy(self.cellregions)
def SimpleGrid_get_bfaces(self):
    return asnumpy(self.bfaces)
def SimpleGrid_get_bfaceregions(self):
    return asnumpy(self.bfaceregions)

shared_ptrSimpleGrid.get_points=SimpleGrid_get_points
shared_ptrSimpleGrid.get_cells=SimpleGrid_get_cells
shared_ptrSimpleGrid.get_cellregions=SimpleGrid_get_cellregions
shared_ptrSimpleGrid.get_bfaces=SimpleGrid_get_bfaces
shared_ptrSimpleGrid.get_bfaceregions=SimpleGrid_get_bfaceregions


}
