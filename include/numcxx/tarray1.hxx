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
        
 

    template<typename T> class TArray1: public TArray<T>
    {
    public:
        using TArray<T>::size;
        using TArray<T>::operator[];
        using TArray<T>::operator=;



        

        // Default constructor.
        TArray1();
        

        /// Construct an empty 1D array.
        ///
        /// \param n0 Size.
        TArray1(index n0);

        /// Construct an 1D array from data pointer
        ///
        /// \param n0 Size.
        /// \param data Pointer to data.
        /// \param deleter Deleter method, \see TArray<T>#_deleter
        TArray1(index n0, T*data, std::function<void(T*p)> deleter);

        /// Construct an 1D array from data pointer
        ///
        /// \param n0 Size.
        /// \param data Pointer to data.
        /// \param datamanager Smartpointer to object managing data.
        /// \see TArray<T>#_datamanager
        TArray1(index n0, T*data, std::shared_ptr<void> datamanager);



        /// Construct 1D Array from std::vector.
        ///
        ///  \param v Smart pointer to vector. A copy (increasin refcount)
        /// is stored as datamanager in created object.
        TArray1(std::shared_ptr<std::vector<T>> v);

        /// Construct 1D Array from std::initializer list.
        TArray1(const std::initializer_list<T> &il );


        /// Construct empty 1D Array
        ///
        /// Mainly for access from python
        static std::shared_ptr<TArray1 <T> > create(index n1);

        /// Construct 1D Array from std::initializer list.
        ///
        /// Static factory method to be used in place of ``std::make_shared``
        /// which has problems to detect the type of the ``initializer_list``.
        static std::shared_ptr<TArray1 <T> > create(const std::initializer_list<T> il);

        /// Create a copy of the array.
        ///
        ///  \return Array of the same size with contents initialized to this
        std::shared_ptr<TArray1 <T> > copy() const;

        /// Create a clone of the array.
        ///
        ///  \return Array of the same size with empty contents.
        std::shared_ptr<TArray1 <T> > clone() const;


        TArray1<T>&  operator=(const TArray1<T> &expr) { assign(*this,expr); return *this;}


        /// Access operator for 1D arrays.
        ///
        /// \param i0  Index of element to be accessed.
        /// \return    Reference to element to be accessed.
        T & operator()(index i0);

        /// C++ style access operator
        ///
        /// \return Plain pointer to start of the array
        T & operator[](index i0);

        /// Element read access.
        /// 
        /// \param i0  Index of element to be accessed.
        /// \return Value of element at i0
        T item(index i0) const;


        /// Element write access.
        /// 
        /// \param i0  Index of element to be accessed.
        /// \param x value to be copied to element at index.
        void itemset(index i0, T x);


        /// Element read access.
        /// 
        /// Getter routine for access from python.
        /// \param i0  Index of element to be accessed.
        /// \return Value of element at i0
        T __getitem__(index i0) const;


        /// Element write access.
        /// 
        /// Setter routine for access from python.
        /// \param i0  Index of element to be accessed.
        /// \param x value to be copied to element at index.
        void __setitem__(index i0,T x);

        /// Print contents of array.
        friend std::ostream & operator<< <T>(std::ostream & s, TArray1<T> &A);

        
    private:
        using TArray<T>::_data;
        using TArray<T>::_idx;
       
    };



}

#include "tarray1-imp.hxx"
#endif
