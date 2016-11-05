#ifndef NUMCXX_TARRAY_H
#define NUMCXX_TARRAY_H

#include <typeinfo>
#include <memory>
#include <stdexcept> 



namespace numcxx
{
    using index= unsigned int;

    /// TArray is the common template base class for arrays and dense matrices
    /// of the numcxx project.
    template<typename T> class TArray
    {
    public:

        typedef T value_type;

        /// Obtain C-pointer of data array.
        /// 
        /// \return Address of C-Array managed by the class which holds
        ///         the data
        T*data() const;

        /// Obtain tensor dimension of array.
        /// 
        /// Tensor dimension is 1 for vectors, 2 for matrices.
        /// \return Dimension.
        index ndim() const;


        /// Obtain size of array.
        /// 
        /// This ist the overall number of elements in the array
        /// \return Size.
        index size() const;

        /// Obtain shape of array for given dimension.
        ///
        /// For 1D arrays, ``shape(0)``  is equivalent to the size
        /// For 2D arrays, ``shape(0)`` is the number of rows
        /// and ``shape(1)`` the number of columns. This corresponds
        /// to the "row major" storage format.
        /// \param dim Tensor dimension.
        /// \return Number of elements in given dimension.
        index shape(const index dim)  const;

        /// Entry value for use in expression templates
        /// 
        const T & operator[](index i0) const;


        /// Add value to all elements.
        ///
        /// \param a Summand for each element.
        void operator+=(const T a);

        /// Subtract value from all elements.
        ///
        /// \param a Value to be subracted from each element.
        void operator-=(const T a);


        /// Multiply all elements by value
        ///
        /// \param a Multiplicator for each element.
        void operator*=(const T a);


        /// Divide each element by value
        ///
        /// \param a Divisor for each element.
        void operator/=(const T a);

        /// Add vector.
        ///
        /// \param a Summand.
        void operator+=(const TArray<T> & a);

        /// Subtract vector.
        ///
        /// \param a Vector to be subtracted
        void operator-=(const TArray<T> & a);


        /// Minimum of array.
        ///
        /// \return Minimum of all elements in array.
        T min() const;

        /// Maximum of array.
        ///
        /// \return Maximum of all elements in array.
        T max() const;

        /// Sum of array.
        ///
        /// \return Sum of all elements in array.
        T sum() const;

        /// Euklidean norm.
        ///
        /// \return Euklidean ( \f$l^2\f$) norm of array.
        T norm2() const;

        /// Sum norm of array.
        ///
        /// \return Sum ( \f$l^1\f$) norm of array.
        T norm1() const;

        /// Maximum norm of array.
        ///
        /// \return Maximum ( \f$l^1\f$) norm of array.
        T normi() const;

        /// Set all elements of array to  given value.
        ///
        /// \param a Assignment value for each element
        void fill(T a);

        /// Fill array with elements form other array.
        ///
        /// \param A Array with given values.
        void fill(const TArray<T> &A);

        /// Fill array with function of elements form other array.
        ///
        /// \param f Function to be performed  on values of \p A.
        /// \param A Array with given values.
        void fill(std::function< T(const T)> f,const TArray<T> & A);

        /// Binary operation on arrays.
        ///
        /// \param f Function performing operation for each index.
        /// \param A First array argument.
        /// \param B Second array argument.
        static void operate(std::function< void ( T& a, T&b)> f, TArray<T> & A, TArray<T> & B);

        /// Ternary operation on arrays.
        ///
        /// \param f Function performing operation for each index.
        /// \param A First array argument.
        /// \param B Second array argument.
        /// \param C Third array argument.
        static void operate(std::function< void ( T& a, T&b,T&c)> f, TArray<T> & A, TArray<T> & B,TArray<T> & C);

        /// Dot product of two arrays.
        ///
        /// \param A First array argument.
        /// \param B Second array argument.
        /// \return  Dot product
        static T dot(const TArray<T>& A,const TArray<T> &B);


    private:

        /// Tensor dimension.
        const index _ndim;

        /// Size of array.
        index _size;
        
        /// Shape vector
        index _shape[3]={0,0,0};

        /// Deleter method.
        /// 
        /// This is the  proper method to be used to  destroy the data
        /// pointer in the array if data manager is null. Depending on
        /// the way  it was constructed,  it may do  nothing, ``free``
        /// the memory, ``delete[]`` the memory, or something else.
        const std::function<void(T*p)> _deleter=nullptr;

        /// Data manager.
        /// 
        /// Smart pointer to some other object managing the data pointer.
        /// If it is not nullptr, the deleter is not called and the 
        /// memory corresponding to the data pointer is freed when the
        /// object behind the data manager is destroyed.
        ///
        /// An example in case is the use of ``shared_ptr<vector> v``
        /// as datamanager and  ``v->data()`` as data pointer
        std::shared_ptr<void>_datamanager =nullptr;

        /// Bounds checker.
        /// 
        /// Bounds check is enabled if the code is compiled with
        /// ``-DNUMCXX_CHECK_BOUNDS``.
        ///
        /// On error it throws ``std::out_of_range`` exception.
        /// \param acc_dim  Dimension to be checked.
        /// \param acc_ndim Tensor dimension to be checked.
        /// \param acc_dim  Index in dimension \p acc_dim to be checked.
        void _check_bounds(index acc_dim, index acc_ndim, index acc_idx) const;


    protected:
        /// Data pointer.
        T* _data=nullptr;
        
        /// Check if all shapes are the same.
        ///
        /// Throw an exception on error
        void _check_square() const;
        
        /// 1D Array index calculation with optional bounds check.
        index _idx(index i0) const;

        /// 2D Array index calculation with optional bounds check.
        index _idx(index i0,index i1)  const;

        /// 3D Array index calculation with optional bounds check.
        index _idx(index i0,index i1,index i2)  const;
    

        /// Construct an empty 1D array.
        ///
        /// \param n0 Size.
        TArray(index n0);

        /// Construct an 1D array from data pointer
        ///
        /// \param n0 Size.
        /// \param data Pointer to data.
        /// \param deleter Deleter method, \see TArray<T>#_deleter
        TArray(index n0, T*data, std::function<void(T*p)> deleter);

        /// Construct an 1D array from data pointer
        ///
        /// \param n0 Size.
        /// \param data Pointer to data.
        /// \param deleter Deleter method.
        /// \see TArray<T>#_datamanager
        TArray(index n0, T*data, std::shared_ptr<void> datamanager);

        /// Construct an empty 2D array.
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        TArray(index n0, index n1);

        /// Construct a 2D array from data pointer
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        /// \param data Pointer to data.
        /// \param deleter Deleter method, \see TArray<T>#_deleter
        TArray(index n0, index n1, T*data,std::function<void(T*p)> deleter);

        /// Construct a 2D array from data pointer
        ///
        /// \param n0 Number of rows
        /// \param n1 Number of columns
        /// \param data Pointer to data.
        /// \param deleter Deleter method.
        /// \see TArray<T>#_datamanager
        TArray(index n0, index n1, T*data,std::shared_ptr<void> datamanager);
        
        /// Default constructor.
        TArray();

        /// Destructor.
        ~TArray();

    };
}


#include "tarray-imp.hxx"
#endif
