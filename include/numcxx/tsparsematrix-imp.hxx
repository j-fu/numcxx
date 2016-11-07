#include <algorithm>
#include <functional>

namespace numcxx
{
    template <typename T>
    inline TSparseMatrix<T>::TSparseMatrix(index n):
        n(n),
        A(std::make_shared<TArray1<T>>(n)),
        JA(std::make_shared<TArray1<int>>(n)),
        IA(std::make_shared<TArray1<int>>(n+1)),
        ext(std::make_shared<Extension>(n))
    {
        auto &ia=*IA;
        auto &ja=*JA;
        auto &a=*A;
        ia=0;
        ja=0;
        a=0.0;
    }

    template <typename T>
    inline void TSparseMatrix<T>::flush()
    {
        if (ext==0) return;
        if (ext->empty()) return;
        int nja_old=JA->size();
        int next=ext->JA->size()-n;
        int nja_new=nja_old+next+ext->next0;
    
        auto new_IA=TArray1<int>::create(n+1);
        auto new_JA=TArray1<int>::create(nja_new);
        auto new_A=TArray1<T>::create(nja_new);
    
        auto &new_ia=*new_IA;
        auto &new_ja=*new_JA;
        auto &new_a=*new_A;
    
        auto &old_ia=*IA;
        auto &old_ja=*JA;
        auto &old_a=*A;
    
        auto &ext_ia=*ext->IA;
        auto &ext_ja=*ext->JA;
        auto &ext_a=*ext->A;
    
        int maxrow_ext=0;
        for (index i=0;i<n;i++)
        {
            int k0=i;
            int lrow=1;
            for (int k=k0; k>=0;k0=k, k=ext_ia[k]) lrow++;
            maxrow_ext=std::max(lrow,maxrow_ext);
        }
    
        std::vector<RowEntry> row(maxrow+maxrow_ext+10);
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
                      [](const RowEntry&e1, const RowEntry &e2)-> bool
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
    
        IA=new_IA;
        JA=new_JA;
        A=new_A;
        maxrow=new_maxrow;
        pattern_has_changed=true;
        ext=0;
        ext=std::make_shared<Extension>(n);
    }


    template <typename T>
    inline void TSparseMatrix<T>::apply(const TArray<T> &u,TArray<T> &v )
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

    template <typename T>
    inline T& TSparseMatrix<T>::operator()(int i, int j) 
    {
        auto &ia=*IA;
        auto &ja=*JA;
        auto &a=*A;
    
        for (int k=ia[i];k<ia[i+1];k++)
            if (ja[k]==j) 
                return a[k];

        return ext->entry(i,j);
    }

    template <typename T>
    inline std::shared_ptr<TMatrix<T>> TSparseMatrix<T>::copy_as_dense()
    {
        flush();
        auto pM=TMatrix<T>::create(n,n);
    
        auto &M=*pM;
        auto &ia=*IA;
        auto &ja=*JA;
        auto &a=*A;
        M=0.0;
        for (int i=0;i<n;i++)
            for (int j=ia[i];j<ia[i+1];j++)
                M(i,ja[j])=a[j];
        return pM;
    }


    template <typename T>
    inline TSparseMatrix<T>::Extension::Extension(index n):
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
    
    template <typename T>
    inline T& TSparseMatrix<T>::Extension::entry(int i, int j)
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
        // std::cout<< i<< " " << j<< std::endl;
        
        if (ja[i]<0)
        {
            ja[i]=j;
            a[i]=0.0;
            next0++;
//                std::cout<< "* "<< i<< std::endl;
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
            
            
    template <typename T>
    inline void TSparseMatrix<T>::Extension::apply(const TArray<T>&u, TArray<T>&v)
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
    
    
    template <typename T>
    inline T TSparseMatrix<T>::Extension::read_entry(int ipart, int i, int j)
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
    
    template <typename T>
    inline bool TSparseMatrix<T>::Extension::empty()
    {
        int next=JA->size()-n;
        if (next+next0==0) return true;
        return false;
    }
    

}
