/////////////////////////////////////////////////
// Implementation of  methods of TArray1<T>

namespace numcxx
{

    
    template<typename T> 
    inline TArray1<T>::TArray1(const std::initializer_list<T> &il ):TArray1(il.size())
    {
        index i=0;
        for (auto x = il.begin() ; x != il.end(); x++,i++) _data[i]= *x;
    }
    
 
    
    
    template<typename T> 
    inline std::ostream & operator << (std::ostream & s, TArray1<T> &A)
    {
        for (index i=0;i<A.size();i++) s <<"[" << i << "]: " <<A(i) << std::endl << std::flush;
        return s;
    }


}
