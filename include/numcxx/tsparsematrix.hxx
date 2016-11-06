#ifndef NUMCXX_TSPARSEMATRIX_H
#define NUMCXX_TSPARSEMATRIX_H



#include <vector>
#include <algorithm>
#include <functional>

#include "tarray1.hxx"

namespace numcxx
{
    template <typename T>
    class TSparseMatrixExtension;

    template<typename T> 
    class TSparseMatrix: TLinOperator<T>
    {
        friend class TSparseMatrixExtension<T>;
        
        std::shared_ptr<TSparseMatrixExtension<T>>  ext;

        /// Row pointers
        std::shared_ptr<IArray1> IA;
        
        /// Column indices
        std::shared_ptr<IArray1> JA;  
        
        /// Entries
        std::shared_ptr<TArray1 <T> > A; 

        bool pattern_has_changed=false;
        const int n;
        int maxrow;
    public:
        TSparseMatrix(index n):
            n(n),
            A(std::make_shared<TArray1<T>>(n)),
            JA(std::make_shared<TArray1<int>>(n)),
            IA(std::make_shared<TArray1<int>>(n+1)),
            ext(std::make_shared<TSparseMatrixExtension<T>>(n))
        {
            auto &ia=*IA;
            auto &ja=*JA;
            auto &a=*A;
            ia=0;
            ja=0;
            a=0.0;
            maxrow=0;
        }
        
        void flush(void)
        {
            // no flush if extension closed
            if (ext==0) return;

            if (!ext->empty())
            {
                ext->flush(*this);
                pattern_has_changed=true;
                ext=0;
                ext=std::make_shared<TSparseMatrixExtension<T>>(n);
            }
        }

        void apply(const TArray1<T> u,TArray1<T> v )
        {
            auto &ia=*IA;
            auto &ja=*JA;
            auto &a=*A;
            for (int i=0;i<n;i++)
            {
                v[i]=0;
                for (int j=ia[i];j<ia[i+1];j++)
                    v[i]+=a[j]*u[ja[j]];
            }
            ext->apply(u,v);
        }
    
        T& operator()(int i, int j) {return entry(i,j);}

        T& entry(int i,int j)
        {
            auto &ia=*IA;
            auto &ja=*JA;
            auto &a=*A;

            for (int k=ia[i];k<ia[i+1];k++)
                if (ja[k]==j) return a[k];

            
            return ext->entry(i,j);
        }

        // std::shared_ptr<TSparseMatrix <T> > copy() const;
        // std::shared_ptr<TSparseMatrix <T> > clone() const;

        // T xentry(const index i, const index j) const {return _data[_idx(i,j)];}


        // template <typename VAL>
        // TSparseMatrix<T>&  operator=(const VAL  &expr)  { assign(*this,expr); return *this;}

        // TSparseMatrix<T>&  operator=(const TSparseMatrix<T> &expr) { assign(*this,expr); return *this;}

        // void apply(const TArray1<T> &u, TArray1<T> &v) const;



    };


    template<typename T> 
    class rowentry
    {
    public: 
        int i;
        T a;
    };
        

    template<typename T> 
    class TSparseMatrixExtension
    {
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

    public:        
        TSparseMatrixExtension(index n):
            n(n),
            A(std::make_shared<std::vector<T>>(n)),
            JA(std::make_shared<std::vector<int>>(n)),
            IA(std::make_shared<std::vector<int>>(n))
        {
            auto &ia=*IA;
            auto &ja=*JA;
            auto &a=*A;

            // initialization for zero main diagonal elements.
            for (int i=0;i<n;i++)
            {
                ia[i]=-1;
                ja[i]=-1;
                a[i]=0;
            }
        }

        T& entry(int i, int j)
        {
            auto &ia=*IA;
            auto &ja=*JA;
            auto &a=*A;

            // Assume "classical" extension scheme
            // without main diagonal entries,
            // and with index shift by -1
            
            // Initial state
            
            // a    0 0 0 0
            // xja  -1 -1 -1 -1
            // xia |-1|-1|-1|-1|
            
    
    
            // Insert 1,1 = 1
            // xja[i]=-1/i triggert Existenz des Elements.
            
            // erstes element der Zeile
            // a    0 1 0 0
            // xja  -1 1 -1 -1
            // xia |-1|-1|-1|-1|
            
            // naechstes Element der Zeile
            
            // a    0 1 0 0 1
            // xja  -1 2 -1 -1 3
            // xia |-1|4|-1|-1|-1


            // do we access row i for the first time ?
            // then we initialize this row
            if (ja[i]==-1)
            {
                ja[i]=j;
                a[i]=0;
                next0++;
                return a[i];
            }
  
    
            // xja /xa entries valid only if xia entry >=0
            int k0=i;
            for (int k=k0; k>=0;k0=k, k=ia[k])
                if (ja[k]==j) // no need to add new entry 
                    // no need to check if xja>=0 as j>=0
                    return a[k]; 
  
            // Access failed  (ia[k0]<0)
            // so we add next entry at end of matrix
            // and put its index into ia[k0]
            
            a.push_back(0);
            ja.push_back(j);
            ia.push_back(0);
            ia[k0]=a.size()-1;
            a[k0]=0.0;
            return a[k0];
        }


        void apply(const TArray1<T>&u, TArray1<T>&v)
        {
            if (empty()) return;
            auto &ia=*IA;
            auto &ja=*JA;
            auto &a=*A;
            
            for (int i=0;i<n; i++)
                for (int k=i; k>=0;k=ia[k])
                {
                    int j=ja[k];
                    if (j>=0)
                        v[i]+=a[k]*u[j];
                }
        }
        

        T read_entry(int ipart, int i, int j)
        {
            auto &ia=*IA;
            auto &ja=*JA;
            auto &a=*A;
            
            
            int k0=i;
            for (int k=k0; k>=0;k0=k, k=ia[k])
                if (ja[k]==j) 
                    return a[k];

            return zero;
        }
        
        bool empty()
        {
            int next=JA->size()-n;
            if (next+next0==0) return true;
            return false;
        }
        
        void  flush(TSparseMatrix<T> & matrix)
        {
            
            if (empty()) return;
            int nja_old=matrix.JA->size();
            int next=JA->size()-n;
            int nja_new=nja_old+next+next0;
            
            auto new_IA=TArray1<int>::create(n+1);
            auto new_JA=TArray1<int>::create(nja_new);
            auto new_A=TArray1<T>::create(nja_new);
            
            auto &new_ia=*new_IA;
            auto &new_ja=*new_JA;
            auto &new_a=*new_A;
            
            auto &old_ia=*matrix.IA;
            auto &old_ja=*matrix.JA;
            auto &old_a=*matrix.A;
            
            auto &ext_ia=*IA;
            auto &ext_ja=*JA;
            auto &ext_a=*A;
            
            int maxrow=0;
            for (index i=0;i<n;i++)
            {
                int k0=i;
                int lrow=1;
                for (int k=k0; k>=0;k0=k, k=ext_ia[k]) lrow++;
                maxrow=std::max(lrow,maxrow);
            }
            
            std::vector<rowentry<T>> row(matrix.maxrow+maxrow+10);
            int new_maxrow=0;
            int j=0;
            int i;
            for(i=0;i<n;i++)
            {
                // put extension entries into row and sort them
                int lxrow=0;
                int k0=i;
                for (int k=k0; k>=0;k0=k, k=ext_ia[k])
                    if (ext_ja[k]>=0)
                    {
                        row[lxrow].i=ext_ja[k];
                        row[lxrow].a=ext_a[k];
                        lxrow++;
                    }

                
                std::sort(row.begin(),row.begin()+lxrow, 
                          [](const rowentry<T>&e1, const rowentry<T> &e2)-> bool
                          {
                              return  (e1.i<e2.i);
                          });
                
                // jointly sort old and new entries into new_ja
                int j0;
                j0=j;
                new_ia[i]=j;
                int irow=0;
                int k=old_ia[i];
                for(;;)
                {
                    if (k<old_ia[i+1] && ((irow>lxrow) ||old_ja[k]<row[irow].i))
                    {
                        new_ja[j]=old_ja[k];
                        new_a[j]=old_a[k];
                        k++;
                        j++;
                        continue;
                    }
                    
                    if (irow<lxrow)
                    {
                        new_ja[j]=row[irow].i;
                        new_a[j]=row[irow].a;
                        irow++;
                        j++;
                        continue;
                    }
                    break;
                }	    
                new_maxrow=std::max(new_maxrow,j-j0);
            }
            new_ia[i]=j;
            
            matrix.IA=new_IA;
            matrix.JA=new_JA;
            matrix.A=new_A;
            matrix.maxrow=new_maxrow;
        }
    };
    
    
    
}


#endif
