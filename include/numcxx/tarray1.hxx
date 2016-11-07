#ifndef NUMCXX_TARRAY1_H
#define NUMCXX_TARRAY1_H

#include <vector>
#include <ostream>
#include <iostream>
#include <type_traits>

#include "tarray.hxx"




namespace numcxx
{

    template<typename T>
    class TArray1;


    template<typename T>
    inline std::ostream & operator << (std::ostream & s, TArray1<T> &A);
        
    
    /// One dimensional array class
    template<typename T> class TArray1: public TArray<T>
    {
    public:
        using TArray<T>::operator[];
        using TArray<T>::size;
        using TArray<T>::operator=;


        // Default constructor.
        TArray1():TArray<T>(){};
        

        /// Construct an empty 1D array.
        ///
        /// \param n0 Size.
        TArray1(index n0):TArray<T>(n0){};

        /// Construct an 1D array from data pointer
        ///
        /// \param n0 Size.
        /// \param data Pointer to data.
        /// \param deleter Deleter method, \see TArray<T>#_deleter
        TArray1(index n0, T*data, std::function<void(T*p)> deleter):TArray<T>(n0,data,deleter){};

        /// Construct an 1D array from data pointer
        ///
        /// \param n0 Size.
        /// \param data Pointer to data.
        /// \param datamanager Smartpointer to object managing data.
        /// \see TArray<T>#_datamanager
        TArray1(index n0, T*data, std::shared_ptr<void> datamanager):TArray<T>(n0,data,datamanager){}; 




        /// Construct 1D Array from std::vector.
        ///
        ///  \param v Smart pointer to vector. A copy (increasin refcount)
        /// is stored as datamanager in created object.
        TArray1(std::shared_ptr<std::vector<T>> v):TArray<T>(v->size(),v->data(),v){};


        /// Construct 1D Array from std::initializer list.
        TArray1(const std::initializer_list<T> &il );

        

        /// Construct empty 1D Array
        ///
        /// Mainly for access from python
        static std::shared_ptr<TArray1 <T> > create(index n1) { return std::make_shared<TArray1 <T> >(n1);}


        /// Construct 1D Array from std::initializer list.
        ///
        /// Static factory method to be used in place of ``std::make_shared``
        /// which has problems to detect the type of the ``initializer_list``.
        static std::shared_ptr<TArray1 <T> > create(const std::initializer_list<T> il){return std::make_shared<TArray1 <T> >(il);}

        /// Create a copy.
        std::shared_ptr<TArray1 <T> > copy() const {return std::make_shared<TArray1 <T>>(*this); }


        /// Create a clone of the array.
        ///
        ///  \return Array of the same size with empty contents.
        std::shared_ptr<TArray1 <T> > clone() const     {  return create(size()); }



        /// Assignment operator
        TArray1<T>&  operator=(const TArray1<T> &expr) {return static_cast<TArray1<T>&>(assign(*this,expr));}

        // Copy constructor
        TArray1(const TArray1<T>& A):TArray1<T>(A.shape(0)){assign(*this,A);}

        /// Access operator for 1D arrays.
        ///
        /// \param i0  Index of element to be accessed.
        /// \return    Reference to element to be accessed.
        T & operator()(index i0)  { return _data[_idx(i0)];};

        /// Element read access.
        /// 
        /// \param i0  Index of element to be accessed.
        /// \return Value of element at i0
        T item(index i0) const  { return _data[_idx(i0)];};


        /// Element write access.
        /// 
        /// \param i0  Index of element to be accessed.
        /// \param x value to be copied to element at index.
        void itemset(index i0, T x) { _data[_idx(i0)]=x;};


        /// Element read access.
        /// 
        /// Getter routine for access from python.
        /// \param i0  Index of element to be accessed.
        /// \return Value of element at i0
        T __getitem__(index i0) const  { return _data[_idx(i0)]; };


        /// Element write access.
        /// 
        /// Setter routine for access from python.
        /// \param i0  Index of element to be accessed.
        /// \param x value to be copied to element at index.
        void __setitem__(index i0,T x)  { _data[_idx(i0)]=x; };

        /// Print contents of array.
        friend std::ostream & operator<< <T>(std::ostream & s, TArray1<T> &A);

        
    private:
        using TArray<T>::_data;
        using TArray<T>::_idx;
       
    };



}

#include "tarray1.ixx"
#endif
