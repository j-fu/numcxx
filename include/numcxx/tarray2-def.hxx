namespace numcxx
{
    
/////////////////////////////////////////////////
// inline methods of TArray2<T>
    
    
    template <typename T> 
    inline TArray2<T>::TArray2():TArray<T>(){};   
    
    template <typename T> 
    inline TArray2<T>::TArray2(index n0, index n1):TArray<T>(n0,n1) {}
    
    template <typename T> 
    inline TArray2<T>::TArray2(index n0, index n1, T*data,std::function<void(T*p)> deleter):
        TArray<T>(n0,n1,data,_deleter(deleter)){}
    
    
    template <typename T> 
    inline TArray2<T>::TArray2(index n0,index n1, T*data, std::shared_ptr<void> datamanager):
        TArray<T>(n0,n1,data,datamanager) {}
    
    template <typename T> 
    inline std::shared_ptr<TArray2 <T> > TArray2<T>::create(index n0,index n1)
    {
        return std::make_shared<TArray2 <T> >(n0,n1);
    }
    
    template <typename T> 
    inline T & TArray2<T>::operator()(index i0, index i1) { return _data[_idx(i0,i1)];};
    
    template <typename T> 
    inline T * TArray2<T>::operator[](index i0) { return &_data[_idx(i0,0)];};

    template <typename T> 
    inline T TArray2<T>::item(index i0,index i1) { return _data[_idx(i0,i1)];};

    template <typename T> 
    inline void TArray2<T>::itemset(index i0, index i1, T x) { _data[_idx(i0,i1)]=x;};

    template <typename T> 
    inline std::shared_ptr<TArray2 <T> > TArray2<T>::copy() const
    {
        auto x=create(shape(0),shape(1));
        x->fill(*this);
        return x;
    }

    template <typename T> 
    inline std::shared_ptr<TArray2 <T> > TArray2<T>::clone() const
    {
        return create(shape(0),shape(1));
    }

    template <typename T> 
    inline std::shared_ptr<TArray1 <T> > const TArray2<T>::__getitem__(index i0) 
    { 
        return std::shared_ptr<TArray1<T>>(new TArray1<T>(shape(1), &_data[_idx(i0,0)], [](T*p){;}));
    } 



    template <typename T> 
    inline TArray2<T>::TArray2(const  std::initializer_list<std::initializer_list<T>> &il ):TArray2(il.size(), il.begin()->size())
    {
        index i=0;
        
        for (auto jl = il.begin() ; jl != il.end(); jl++,i++)
        {
            index j=0;
            for (auto x = jl->begin() ; x != jl->end(); x++,j++) 
                _data[_idx(i,j)]= *x;
        }
    }
    

    template <typename T> 
    inline std::shared_ptr<TArray2 <T> > TArray2<T>::create(const  std::initializer_list<std::initializer_list<T>> &il)
    {
        return std::make_shared<TArray2<T>> (il);
    }
        

    
    
    template <typename T> 
    inline std::ostream & operator << (std::ostream & s, TArray2<T> &A)
    {
        s << "    ";
        for (index j=0;j<A.shape(1);j++) 
            s << "[" << j << "]     ";
        s<< std::endl;
        for (index i=0;i<A.shape(0);i++) 
        {
            s << "[" << i << "]: ";
            for (index j=0;j<A.shape(1);j++) 
                s << A(i,j) << "   ";
            s<< std::endl;
            }
        s << std::flush;
        return s;
    }
    
    template <typename T> 
    inline std::ostream & operator << (std::ostream & s, const std::shared_ptr<TArray2<T>>&A)
    {
        return operator<<(s,*A);
    }

}
