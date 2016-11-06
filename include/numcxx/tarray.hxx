#ifndef NUMCXX_TARRAY_H
#define NUMCXX_TARRAY_H

#include <typeinfo>
#include <memory>
#include <stdexcept> 



namespace numcxx
{
    using index= unsigned int;

    class ExpressionBase
    {
    };

    /// TArray is the common template base class for arrays and dense matrices
    /// of the numcxx project.
    template<typename T> class TArray: public ExpressionBase
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

        /// Const reference to entry for use in expression templates
        const T & operator[](index i0) const;
        
        template <typename VAL>
        TArray<T>&  operator=(const VAL  &expr)  { assign(*this,expr); return *this;}


        /// Add value to all elements.
        ///
        /// \param a Summand for each element.
        template <typename VAL>
        void operator+=(const VAL & a) {xadd(*this,a);}

        /// Subtract value from all elements.
        ///
        /// \param a Value to be subracted from each element.
        template <typename VAL>
        void operator-=(const VAL& a) {xsub(*this,a);}


        /// Multiply all elements by value
        ///
        /// \param a Multiplicator for each element.
        template <typename VAL>
        void operator*=(const VAL& a) {xmul(*this,a);}


        /// Divide each element by value
        ///
        /// \param a Divisor for each element.
        template <typename VAL>
        void operator/=(const VAL & a) {xdiv(*this,a);}



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

    /// Maximum norm of array.
    ///
    /// \return Maximum ( \f$l^\infty\f$) norm.
    template <typename A> double normi(const A& a);

    /// Sum norm.
    ///
    /// \return Sum ( \f$l^1\f$) norm.
    template <typename A> double norm1(const A& a);

    /// Euklidean norm.
    ///
    /// \return Euklidean ( \f$l^2\f$) norm 
    template <typename A> double norm2(const A& a);

    /// Dot product.
    ///
    /// \param A First argument.
    /// \param B Second argument.
    /// \return  Dot product
    template <typename A, typename B> double dot(const A& a, const B&b);

    
    /// Minimum of array.
    ///
    /// \return Minimum of all elements in array.
    template <typename A> double min(const A&a);
    
    /// Maximum of array.
    ///
    /// \return Maximum of all elements in array.
    template <typename A>  double max(const A&a);
    
    /// Sum of array.
    ///
    /// \return Sum of all elements in array.
    template <typename A>  double sum(const A&a);
    
    template<typename T> class TLinSolver
    {
    public:
        TLinSolver(){};
        virtual void solve( TArray<T> & sol,  const TArray<T> & rhs) const {};
        virtual void update(void){};
    };

    template<typename T> class TLinOperator
    {
    public:
        TLinOperator(){};
        virtual void apply( const TArray<T> & sol,   TArray<T> & rhs) const {};
    };

    
    
}


#include "tarray-imp.hxx"
#endif
