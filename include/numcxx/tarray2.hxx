#ifndef NUMCXX_TARRAY2_H
#define NUMCXX_TARRAY2_H

#include <ostream>
#include "tarray.hxx"


namespace  numcxx
{
    template<typename T>
    class TArray2;

    template<typename T>
    inline std::ostream & operator << (std::ostream & s, TArray2<T> &A);

    template<typename T> 
    class TArray2: public TArray<T>
    {
    public:
        using TArray<T>::size;
        using TArray<T>::shape;
        using TArray<T>::operator[];
        using TArray<T>::operator=;
        
        // Default constructor.
        TArray2();

        /// Construct an empty 2D array.
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        TArray2(index n0, index n1);

        /// Construct a 2D array from data pointer
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        /// \param data Pointer to data.
        /// \param deleter Deleter method, \see TArray<T>#_deleter
        TArray2(index n0, index n1, T*data,std::function<void(T*p)> deleter);

        /// Construct a 2D array from data pointer
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        /// \param data Pointer to data.
        /// \param deleter Deleter method.
        /// \see TArray<T>#_datamanager
        TArray2(index n0, index n1, T*data,std::shared_ptr<void> datamanager);

        /// Construct 2D Array from std::initializer list.
        TArray2(const  std::initializer_list<std::initializer_list<T>> &il );

        TArray2<T>&  operator=(const TArray2<T> &expr) { assign(*this,expr); return *this;}


        /// Construct empty 2D Array
        ///
        /// Mainly for access from python
        static std::shared_ptr<TArray2 <T> > create(index n0,index n1);

        /// Construct 2D Array from std::initializer list.
        static std::shared_ptr<TArray2 <T> > create(const  std::initializer_list<std::initializer_list<T>> &il);

        /// Create a copy of the array.
        ///
        ///  \return Array of the same size with contents initialized to this
        std::shared_ptr<TArray2 <T> > copy() const;

        /// Create a clone of the array.
        ///
        ///  \return Array of the same size with empty contents.
        std::shared_ptr<TArray2 <T> > clone() const;


        /// Access operator for 2D arrays.
        ///
        /// \param i0  Row index of element to be accessed.
        /// \param i0  Column index of element to be accessed.
        /// \return    Reference to element to be accessed.
        T & operator()(index i0, index i1);


        /// C++ style access operator
        /// 
        /// \return Plain pointer to i0-th  row of the array
        T * operator[](index i0);


        /// Element read access.
        /// 
        /// \param i0  Row index of element to be accessed.
        /// \param i1  Column index of element to be accessed.
        /// \return Value of element at (i0,i1)
        T item(index i0,index i1);


        /// Element write access.
        /// 
        /// \param i0  Row index of element to be accessed.
        /// \param i1  Column index of element to be accessed.
        /// \param x value to be copied to element at index.
        void itemset(index i0, index i1, T x);

        /// Getter routine for access from python.
        /// 
        /// This access is rather expensive, as it constructs
        /// a smart pointer to the row. 
        /// \param i0  row index of element to be accessed.
        /// \return Smart pointer to i0-th row. 
        std::shared_ptr<TArray1 <T> > const __getitem__(index i0);

        /// Print contents of array.
        friend std::ostream & operator<< <T>(std::ostream & s, TArray2<T> &A);

    protected:
        using TArray<T>::_data;
        using TArray<T>::_idx;
    };

}
#include "tarray2-imp.hxx"

#endif
