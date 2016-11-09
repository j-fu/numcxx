#include <algorithm>
#include <functional>

namespace numcxx
{


    template <typename T>
    class TSparseMatrix<T>::RowEntry
    {
    public: 
        index i;
        T a;
    };
    
    template <typename T>
    class TSparseMatrix<T>::Extension
    {
    public:        
        /// number of matrix rows
        const index n;               
        
        
        /// rowpointer list 
        std::shared_ptr<std::vector<int>> pIA;     
        
        /// column indices 
        std::shared_ptr<std::vector<int>> pJA;     
        
        /// values
        std::shared_ptr<std::vector <T> > pA; 
        
        /// zero value
        const T zero=0;        
        
        /// Constructor
        Extension(index n);
        
        /// Return refrence to  entry, if not existent, create
        T& entry(int i, int j);
        
        void apply(const TArray<T>&U, TArray<T>&V);
        
        int next;
        
        T read_entry(int i, int j);
        
        bool empty();
    };
    





    template <typename T>
    inline TSparseMatrix<T>::TSparseMatrix(index n):
        n(n),
        pA(std::make_shared<TArray1<T>>(n)),
        pJA(std::make_shared<TArray1<int>>(n)),
        pIA(std::make_shared<TArray1<int>>(n+1))
    {
        auto &IA=*pIA;
        auto &JA=*pJA;
        auto &A=*pA;
        IA=0;
        JA=0;
        A=0.0;
    }
    
    template <typename T>
    inline TSparseMatrix<T>::TSparseMatrix(const  std::initializer_list<std::initializer_list<T>> &il):
        TSparseMatrix(il.size())
    {
        index i=0;
        for (auto jl = il.begin() ; jl != il.end(); jl++,i++)
        {
            index j=0;
            for (auto x = jl->begin() ; x != jl->end(); x++,j++) 
                this->operator()(i,j)= *x;

        }
        flush();
    }

    template <typename T>
    inline std::shared_ptr<TSparseMatrix <T> > TSparseMatrix<T>::create(const  std::initializer_list<std::initializer_list<T>> &il)
    {
        return std::make_shared<TSparseMatrix <T>>(il);
    }

    template <typename T>
    inline void TSparseMatrix<T>::flush()
    {
        if (pExt==nullptr) return;
        if (pExt->empty()) return;
        int nJA_old=0;
        if (!empty())
            nJA_old=pJA->size();
        int next=pExt->next;
        int nJA_new=nJA_old+next;
    
        auto pNew_IA=TArray1<int>::create(n+1);
        auto pNew_JA=TArray1<int>::create(nJA_new);
        auto pNew_A=TArray1<T>::create(nJA_new);
    
        auto &New_IA=*pNew_IA;
        auto &New_JA=*pNew_JA;
        auto &New_A=*pNew_A;
    
        auto &Old_IA=*pIA;
        auto &Old_JA=*pJA;
        auto &Old_A=*pA;
    
        auto &Ext_IA=*pExt->pIA;
        auto &Ext_JA=*pExt->pJA;
        auto &Ext_A=*pExt->pA;
    
        int maxrow_ext=0;
        for (index i=0;i<n;i++)
        {
            int k0=i;
            int lrow=1;
            for (int k=k0; k>=0;k0=k, k=Ext_IA[k]) lrow++;
            maxrow_ext=std::max(lrow,maxrow_ext);
        }
    
        std::vector<RowEntry> row(maxrow+maxrow_ext+10);
        int New_maxrow=0;
        int j=0;
        int i;
        for(i=0;i<n;i++)
        {
            // put extension entries into row and sort them
            int lxrow=0;
            int k0=i;
            for (int k=k0; k>=0;k0=k, k=Ext_IA[k])
                if (Ext_JA[k]>=0)
                {
                    row[lxrow].i=Ext_JA[k];
                    row[lxrow].a=Ext_A[k];
                    lxrow++;
                }
        
        
            std::sort(row.begin(),row.begin()+lxrow, 
                      [](const RowEntry&e1, const RowEntry &e2)-> bool
                      {
                          return  (e1.i<e2.i);
                      });
        
            // jointly sort Old and New entries into New_JA
            int j0;
            j0=j;
            New_IA[i]=j;
            int irow=0;
            int k=Old_IA[i];
        
            for(;;)
            {
                if (k<Old_IA[i+1] && ((irow>lxrow) ||Old_JA[k]<row[irow].i))
                {
                    New_JA[j]=Old_JA[k];
                    New_A[j]=Old_A[k];
                    k++;
                    j++;
                    continue;
                }
            
                if (irow<lxrow)
                {
                    New_JA[j]=row[irow].i;
                    New_A[j]=row[irow].a;
                    irow++;
                    j++;
                    continue;
                }
                break;
            }	    
            New_maxrow=std::max(New_maxrow,j-j0);
        }
        New_IA[i]=j;
    
        pIA=pNew_IA;
        pJA=pNew_JA;
        pA=pNew_A;
        maxrow=New_maxrow;
        _pattern_changed=true;
        _first_flush_done=true;
        pExt=nullptr;
    }


    template <typename T>
    inline void TSparseMatrix<T>::apply(const TArray<T> &u,TArray<T> &v ) const
    {
        v.resize(u.size());
        auto &IA=*pIA;
        auto &JA=*pJA;
        auto &A=*pA;
        for (int i=0;i<n;i++)
        {
            v[i]=0;
            for (int j=IA[i];j<IA[i+1];j++)
                v[i]+=A[j]*u[JA[j]];
        }
        if (pExt!=nullptr)
            pExt->apply(u,v);
    }

    template <typename T>
    inline T& TSparseMatrix<T>::operator()(int i, int j) 
    {
        auto &IA=*pIA;
        auto &JA=*pJA;
        auto &A=*pA;
    
        for (int k=IA[i];k<IA[i+1];k++)
            if (JA[k]==j) 
                return A[k];
        
        if (pExt==nullptr) pExt=std::make_shared<Extension>(n);
        return pExt->entry(i,j);
    }

    template <typename T>
    inline std::shared_ptr<TMatrix<T>> TSparseMatrix<T>::copy_as_dense()
    {
        flush();
        auto pM=TMatrix<T>::create(n);
    
        auto &M=*pM;
        auto &IA=*pIA;
        auto &JA=*pJA;
        auto &A=*pA;
        M=0.0;
        for (int i=0;i<n;i++)
            for (int j=IA[i];j<IA[i+1];j++)
                M(i,JA[j])=A[j];
        return pM;
    }


        /// Row entry for flush method: we need to sort
        /// both i and a values at the same time.



    template <typename T>
    inline TSparseMatrix<T>::Extension::Extension(index n):
        n(n),
        pA(std::make_shared<std::vector<T>>(n)),
        pJA(std::make_shared<std::vector<int>>(n)),
        pIA(std::make_shared<std::vector<int>>(n))
    {
        auto &IA=*pIA;
        auto &JA=*pJA;
        auto &A=*pA;
        
        // initIAlization for zero main dIAgonal elements.
        for (int i=0;i<n;i++)
        {
            IA[i]=-1;
            JA[i]=-1;
            A[i]=0;
        }
    }
    
    template <typename T>
    inline T& TSparseMatrix<T>::Extension::entry(int i, int j)
    {
        auto &IA=*pIA;
        auto &JA=*pJA;
        auto &A=*pA;
        
        // Assume "classical" extension scheme
        // without main dIAgonal entries,
        // and with index shift by -1
        
        // InitIAl state
        
        // a    0 0 0 0
        // xJA  -1 -1 -1 -1
        // xIA |-1|-1|-1|-1|
        
        
        
        // Insert 1,1 = 1
        // xJA[i]=-1/i triggert Existenz des Elements.
        
        // erstes element der Zeile
        // a    0 1 0 0
        // xJA  -1 1 -1 -1
        // xIA |-1|-1|-1|-1|
        
        // naechstes Element der Zeile
        
        // a    0 1 0 0 1
        // xJA  -1 2 -1 -1 3
        // xIA |-1|4|-1|-1|-1
        
        
        // do we access row i for the first time ?
        // then we initIAlize this row
        
        if (JA[i]<0)
        {
            JA[i]=j;
            A[i]=0.0;
            next++;
            return A[i];
        }
        
                
        // xJA /xa entries valid only if xIA entry >=0
        int k0=i;
        for (int k=k0; k>=0;k0=k, k=IA[k])
        {

            if (JA[k]==j) // no need to add new entry 
                return A[k]; 
        }
        // Access failed  (IA[k0]<0)
        // so we add next entry at end of matrix
        // and put its index into IA[k0]
        
        next++;
        A.push_back(0);
        JA.push_back(j);
        IA.push_back(-1);
        IA[k0]=A.size()-1;
        return A[IA[k0]];
    }
            
            
    template <typename T>
    inline void TSparseMatrix<T>::Extension::apply(const TArray<T>&u, TArray<T>&v)
    {
        
        if (empty()) return;
        auto &IA=*pIA;
        auto &JA=*pJA;
        auto &A=*pA;
        
        for (int i=0;i<n; i++)
            for (int k=i; k>=0;k=IA[k])
            {
                int j=JA[k];
                if (j>=0)
                    v[i]+=A[k]*u[j];
            }
    }
    
    
    template <typename T>
    inline T TSparseMatrix<T>::Extension::read_entry(int i, int j)
    {
        auto &IA=*pIA;
        auto &JA=*pJA;
        auto &A=*pA;
        
        
        int k0=i;
        for (int k=k0; k>=0;k0=k, k=IA[k])
            if (JA[k]==j) 
                return A[k];
        
        return zero;
    }
    
    template <typename T>
    inline bool TSparseMatrix<T>::Extension::empty()
    {
        if (next==0) return true;
        return false;
    }
    

}
