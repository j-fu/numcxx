/////////////////////////////////////////////////
// Implementation of  methods of TArray1<T>

namespace numcxx
{


    template <typename T> 
    inline TArray1<T>::TArray1():TArray<T>(){};   

    template <typename T> 
    inline TArray1<T>::TArray1(index n0):TArray<T>(n0){}

    template <typename T> 
    inline TArray1<T>::TArray1(index n0, T*data,std::function<void(T*p)> deleter):
        TArray<T>(n0,data,deleter){}
    
    template <typename T> 
    inline TArray1<T>::TArray1(index n0, T*data, std::shared_ptr<void> datamanager):
        TArray<T>(n0,data,datamanager){}; 

    template <typename T> 
    inline TArray1<T>::TArray1(std::shared_ptr<std::vector<T>> v):
        TArray<T>(v->size(),v->data(),v){};


    

    template <typename T> 
    inline std::shared_ptr<TArray1 <T> > TArray1<T>::create(index n1)
    {
        return std::make_shared<TArray1 <T> >(n1);
    }
    
    template<typename T> 
    inline TArray1<T>::TArray1(const std::initializer_list<T> &il ):TArray1(il.size())
    {
        index i=0;
        for (auto x = il.begin() ; x != il.end(); x++,i++) _data[i]= *x;
    }
    
    template<typename T>
    inline std::shared_ptr<TArray1 <T> > TArray1<T>::create(const std::initializer_list<T> il)
    {
        return std::make_shared<TArray1 <T> >(il);
    }
    
    template <typename T> 
    inline T & TArray1<T>::operator()(index i0) { return _data[_idx(i0)];};

    template <typename T> 
    inline T & TArray1<T>::operator[](index i0) { return _data[_idx(i0)];};
    
    template <typename T> 
    inline T TArray1<T>::item(index i0) const { return _data[_idx(i0)];};

    template <typename T> 
    inline void TArray1<T>::itemset(index i0, T x) { _data[_idx(i0)]=x;};

    template <typename T> 
    inline std::shared_ptr<TArray1 <T> > TArray1<T>::copy() const
    {
        auto x=create(size());
        *x=*this;
        return x;
    }

    template <typename T> 
    inline std::shared_ptr<TArray1 <T> > TArray1<T>::clone() const
    {
        return create(size());
    }

    template <typename T> 
    inline T TArray1<T>::__getitem__(index i) const { return _data[_idx(i)]; } 

    template <typename T> 
    inline void TArray1<T>::__setitem__(index i,T d) { _data[_idx(i)]=d; } 

    template<typename T> 
    inline std::ostream & operator << (std::ostream & s, const std::shared_ptr<TArray1<T>>&A)
    {
        return operator<<(s,*A);
    }
    
    
    template<typename T> 
    inline std::ostream & operator << (std::ostream & s, TArray1<T> &A)
    {
        for (index i=0;i<A.size();i++) s <<"[" << i << "]: " <<A(i) << std::endl << std::flush;
        return s;
    }


}
