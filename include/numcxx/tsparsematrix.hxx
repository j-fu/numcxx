#ifndef NUMCXX_TSPARSEMATRIX_H
#define NUMCXX_TSPARSEMATRIX_H

#include <vector>
#include "tarray1.hxx"
#include "tmatrix.hxx"

namespace numcxx
{
    
    template<typename T> 
    class TSolverUMFPACK;

    /// Sparse matrix class using CRS storage scheme
    template<typename T> 
    class TSparseMatrix: TLinOperator<T>
    {
        friend class TSolverUMFPACK<T>;
    public:
        TSparseMatrix(index n);
        TSparseMatrix(const  std::initializer_list<std::initializer_list<T>> &il);
        void  flush();
        void apply(const TArray<T> &U,TArray<T> &V );
        T& operator()(int i, int j);
        std::shared_ptr<TMatrix<T>> copy_as_dense();
        static std::shared_ptr<TSparseMatrix <T> > create(const  std::initializer_list<std::initializer_list<T>> &il);


        index shape(int idim) {return n;}
        bool pattern_changed(){return _pattern_changed;}
        void pattern_changed(bool chg) {_pattern_changed=chg;};

        // std::shared_ptr<TSparseMatrix <T> > copy() const;
        // std::shared_ptr<TSparseMatrix <T> > clone() const;
        // T xentry(const index i, const index j) const {return _data[_idx(i,j)];}
        // template <typename VAL>
        // TSparseMatrix<T>&  operator=(const VAL  &expr)  { assign(*this,expr); return *this;}
        // TSparseMatrix<T>&  operator=(const TSparseMatrix<T> &expr) { assign(*this,expr); return *this;}
        // void apply(const TArray1<T> &u, TArray1<T> &v) const;

    private:
        class RowEntry
        {
        public: 
            int i;
            T a;
        };
        
        /// Row pointers
        std::shared_ptr<IArray1> pIA;
        
        /// Column indices
        std::shared_ptr<IArray1> pJA;  
        
        /// Entries
        std::shared_ptr<TArray1 <T> > pA; 
        
        bool _pattern_changed=false;
        const int n;
        int maxrow=0;

        class Extension
        {
        public:        
            /// number of matrix rows
            const index n;               
            
            int next0=0;
            
            /// rowpointer list 
            std::shared_ptr<std::vector<int>> pIA;     
            
            /// column indices 
            std::shared_ptr<std::vector<int>> pJA;     
            
            /// values
            std::shared_ptr<std::vector <T> > pA; 
            
            /// zero value
            const T zero=0;        
            
            Extension(index n);
            
            T& entry(int i, int j);
            void apply(const TArray<T>&U, TArray<T>&V);
            
            
            T read_entry(int i, int j);
            
            bool empty();
        };

        std::shared_ptr<Extension>  pExt;
        
    };


        

    
    
    
}

#include "tsparsematrix-imp.hxx"
#endif
