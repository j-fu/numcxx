/////////////////////////////////////////////////
// inline methods of TArray2<T>

namespace numcxx
{


    template <typename T> 
    inline TArray2<T>::TArray2(const  std::initializer_list<std::initializer_list<T>> &il):
        TArray2(il.size(), il.begin()->size())
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
    

}
