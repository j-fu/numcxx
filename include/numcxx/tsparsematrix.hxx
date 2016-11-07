#ifndef NUMCXX_TSPARSEMATRIX_H
#define NUMCXX_TSPARSEMATRIX_H

#include <vector>
#include "tarray1.hxx"
#include "tmatrix.hxx"

namespace numcxx
{
    
    /// Sparse matrix class using CRS storage scheme
    template<typename T> 
    class TSparseMatrix: TLinOperator<T>
    {
    public:
        TSparseMatrix(index n);
        void  flush();
        void apply(const TArray<T> &u,TArray<T> &v );
        T& operator()(int i, int j);
        std::shared_ptr<TMatrix<T>> copy_as_dense();



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
        std::shared_ptr<IArray1> IA;
        
        /// Column indices
        std::shared_ptr<IArray1> JA;  
        
        /// Entries
        std::shared_ptr<TArray1 <T> > A; 
        
        bool pattern_has_changed=false;
        const int n;
        int maxrow=0;

        class Extension
        {
        public:        
            /// number of matrix rows
            const index n;               
            
            int next0=0;
            
            /// rowpointer list 
            std::shared_ptr<std::vector<int>> IA;     
            
            /// column indices 
            std::shared_ptr<std::vector<int>> JA;     
            
            /// values
            std::shared_ptr<std::vector <T> > A; 
            
            /// zero value
            const T zero=0;        
            
            Extension(index n);
            
            T& entry(int i, int j);
            void apply(const TArray<T>&u, TArray<T>&v);
            
            
            T read_entry(int ipart, int i, int j);
            
            bool empty();
        };

        std::shared_ptr<Extension>  ext;
        
    };


        

    
    
    
}

#include "tsparsematrix-imp.hxx"
#endif
