#ifndef NUMCXX_TARRAY2_H
#define NUMCXX_TARRAY2_H

#include "tarray.hxx"


namespace  numcxx
{
    /// Two-dimensional array class
    template<typename T> 
    class TArray2: public TArray<T>, public ExpressionBase
    {
    public:
        using TArray<T>::size;
        using TArray<T>::shape;
        using TArray<T>::operator[];
        using TArray<T>::operator=;
        
        // Construct zero size array.
        TArray2():TArray<T>(){};;

        /// Construct an empty 2D array.
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        TArray2(index n0, index n1):TArray<T>(n0,n1) {};

        /// Construct a 2D array from data pointer
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        /// \param data Pointer to data.
        /// \param deleter Deleter method, \see TArray<T>#_deleter
        TArray2(index n0, index n1, T*data,std::function<void(T*p)> deleter):TArray<T>(n0,n1,data,_deleter(deleter)){};


        /// Construct a 2D array from data pointer
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        /// \param data Pointer to data.
        /// \param deleter Deleter method.
        /// \see TArray<T>#_datamanager
        TArray2(index n0, index n1, T*data,std::shared_ptr<void> datamanager):TArray<T>(n0,n1,data,datamanager) {};


        /// Construct 2D Array from std::initializer list.
        TArray2(const  std::initializer_list<std::initializer_list<T>> &il ):TArray<T>(il){};


        /// Copy constructor
        TArray2(const TArray2<T>& A):TArray2<T>(A.shape(0),A.shape(1)){assign(*this,A);}


        /// Assignment operator
        TArray2<T>&  operator=(const TArray2<T> &expr) { return static_cast<TArray2<T>&>(assign(*this,expr));}



        /// Construct empty 2D Array
        ///
        /// Mainly for access from python
        static std::shared_ptr<TArray2 <T> > create(index n0,index n1) { return std::make_shared<TArray2 <T> >(n0,n1);}


        /// Construct 2D Array from std::initializer list.
        static std::shared_ptr<TArray2 <T> > create(const  std::initializer_list<std::initializer_list<T>> &il) { return std::make_shared<TArray2 <T> >(il);}

        /// Create a copy of the array.
        ///
        ///  \return Array of the same size with contents initialized to this
        std::shared_ptr<TArray2 <T> > copy() const { return std::make_shared<TArray2 <T> >(*this);}

        /// Create a clone of the array.
        ///
        ///  \return Array of the same size with empty contents.
        std::shared_ptr<TArray2 <T> > clone() const  { return create(shape(0),shape(1));}



        /// Element read access.
        /// 
        /// \param i0  Row index of element to be accessed.
        /// \param i1  Column index of element to be accessed.
        /// \return Value of element at (i0,i1)
        T item(index i0,index i1)  { return _data[_idx(i0,i1)];};


        /// Element write access.
        /// 
        /// \param i0  Row index of element to be accessed.
        /// \param i1  Column index of element to be accessed.
        /// \param x value to be copied to element at index.
        void itemset(index i0, index i1, T x)  { _data[_idx(i0,i1)]=x;};

        /// Getter routine for access from python.
        /// 
        /// This access is rather expensive, as it constructs
        /// a smart pointer to the row. 
        /// \param i0  row index of element to be accessed.
        /// \return Smart pointer to i0-th row. 
        std::shared_ptr<TArray1 <T> > const __getitem__(index i0){ return std::shared_ptr<TArray1<T>>(new TArray1<T>(shape(1), &_data[_idx(i0,0)], [](T*p){;}));}


        bool is_matrix(){return false;}

    protected:
        using TArray<T>::_data;
        using TArray<T>::_idx;
    };

}

#endif
