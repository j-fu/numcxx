/* -*- C -*-  (not really, but good for syntax highlighting) */
%module numcxx

%{
#define SWIG_FILE_WITH_INIT
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numcxx.h"
  namespace numcxx
  {
    class NumpyProxy: public TArrayBase
    {
      PyObject *o;
    public:
    NumpyProxy(PyObject *o): o(o) { Py_INCREF(o);};
      ~NumpyProxy() {Py_DECREF(o);};
    };

    std::shared_ptr<DArray1 > NumpyAsDArray1(PyObject *o,double* v, int n1)
      {
          auto proxy=std::shared_ptr<TArrayBase>(new NumpyProxy(o));
          return std::shared_ptr<DArray1 >(new DArray1 (n1,v,proxy));
      }
    std::shared_ptr<DArray2 > NumpyAsDArray2(PyObject *o,double* v, int n1, int n2)
      {
          return std::shared_ptr<DArray2>(new DArray2 (n1, n2, v,std::shared_ptr<TArrayBase>(new NumpyProxy(o))));
      }

    void DArray1AsNumpy(std::shared_ptr<DArray1 > a,PyObject *o,double** ARGOUTVIEW_ARRAY1, int* DIM1)
    {
        *DIM1=a->size;
        *ARGOUTVIEW_ARRAY1=a->rawdata();
    }
    void DArray2AsNumpy(std::shared_ptr<DArray2 > a,PyObject *o,double** ARGOUTVIEW_ARRAY2, int* DIM1,int* DIM2)
    {
        *DIM1=a->shape[0];
        *DIM2=a->shape[1];
        *ARGOUTVIEW_ARRAY2=a->rawdata();
    }
  }

%}



%include <std_vector.i>
%include <numpy.i>
%clear(double** ARGOUTVIEW_ARRAY1, int* DIM1);
%clear(int** ARGOUTVIEW_ARRAY1, int* DIM1);





/* Typemaps for pdelib <-> numpy array communication.
   
   This typemap  file needs to be included after the original,
   unmodified numpy.i.

   The general idea is to catch  the python objects in order to 
   use Py_INCREF/Py_DECREF in the proper places.
   
*/

%define %pdelib_numpy_typemaps(DATA_TYPE, DATA_TYPECODE, DIM_TYPE)

 /***************************************************************/
 /***************************************************************/
 /***************************************************************/
 /*  
     Numpy array as pdelib array.
 */

 /***************************************************************
   Typemap suite for (PyObject *O, DATA_TYPE* INPLACE_ARRAY1, int DIM1)
   for numpy arrays as pdelib arrays
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
   for numpy arrays as pdelib arrays
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
    pdelib array as  numpy: 
    
    Concerning PyArray_New  instead of PyArray_SimpleNewFromData, see the source on
    https://github.com/numpy/numpy/blob/master/numpy/core/src/multiarray/ctors.c:

    The obj (last) parameter is increfed and stored, so we pass the python proxy
    of the pdelib array.
*/

/***************************************************************
   Typemap suite for (PyObject *O, DATA_TYPE** ARGOUTVIEW_ARRAY1, int* DIM1)
   for pdelib array as  numpy.
   
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
   for pdelib array as  numpy.
   
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

%enddef    /* %pdelib_numpy_typemaps() macro */

%pdelib_numpy_typemaps(double, NPY_DOUBLE, int)
%pdelib_numpy_typemaps(int,NPY_INT, int)



%init %{
  import_array();
%}


/* template <class A> class std::shared_ptr */
/* { */
/*     friend A; */
/* //    aptr(){}; */
/* }; */

template<class A> struct std::shared_ptr
{
    A* operator->() const;
};
  

namespace numcxx
{
    typedef unsigned int index;


    class DArray1
    {
    public:
        const index ndim;
        index size;
        void fill(double x);
        std::vector<index> shape;
        
        static std::shared_ptr<DArray1> create(index n0);
        double item(index i0);
        void itemset(index i0, double x);
        
        std::shared_ptr<DArray1> copy();
        double __getitem__(index i) const;
        void __setitem__(index i,double);
    };
    
    class DArray2
    {
    public:
        const index ndim;
        index size;
        void fill(double x);
        std::vector<index> shape;
        
        static std::shared_ptr<DArray2> create(index n0, index n1);
        double item(index i0,index i1);
        void itemset(index i0,index i1, double x);
        
        std::shared_ptr<DArray2> copy();
        std::shared_ptr<DArray1>   __getitem__(int i);
    };
    
    
    std::shared_ptr<DArray1> NumpyAsDArray1(PyObject *O, double* INPLACE_ARRAY1, int DIM1);
    std::shared_ptr<DArray2> NumpyAsDArray2(PyObject *O, double* INPLACE_ARRAY2, int DIM1, int DIM2);
    
    void DArray1AsNumpy(std::shared_ptr<DArray1> a,PyObject *O,double** ARGOUTVIEW_ARRAY1, int* DIM1);
    void DArray2AsNumpy(std::shared_ptr<DArray2> a,PyObject *O,double** ARGOUTVIEW_ARRAY2, int* DIM1,  int* DIM2);
    
    
}

%template(shared_ptrDArray1) std::shared_ptr<numcxx::DArray1 >;
%template(shared_ptrDArray2) std::shared_ptr<numcxx::DArray2 >;
// for shape etc.
%template(vectoru) std::vector<unsigned int>;

%pythoncode{
from  numpy import ndarray,array,float64
def asnumcxx(proto):
    if type(proto)==ndarray:
        a=proto
    elif type(proto)==list:
        a=array(proto)
    else:
        raise TypeError("asnumcxx: wrong type of input -- must be numpy.ndarray or list.")

    if a.dtype==float64 and a.ndim==1:
        return NumpyAsDArray1(a)          
    elif a.dtype==float64 and a.ndim==2:
        return NumpyAsDArray2(a)           

    raise TypeError("asnumcxx: error")


def asnumpy(a):
    if a.ndim==1:
        return DArray1AsNumpy(a)
    elif a.ndim==2:
        return DArray2AsNumpy(a)

    raise RuntimeError("asnumpy: error")


 }
