#ifndef NUMCXX_EXPRESSION_HXX
#define NUMCXX_EXPRESSION_HXX

namespace numcxx
{
    ////////////////////////////////////////////////////////////////////////////////////
    /// Expression template for   A+B
    template<typename A, typename B,
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    class AdditionExpression
    {
        const A& a;
        const B& b;

    public:

        typedef typename A::value_type value_type;

        AdditionExpression(const A& a, const B& b):a(a), b(b){}
        
        unsigned int size() const { return a.size(); }
            
        const value_type operator[](const unsigned int i) const
        { 
            return a[i] + b[i];
        }
    };
        
    template<typename A, typename B, 
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    const AdditionExpression<A,B> operator+(const A& a, const B& b)
    {
        return AdditionExpression<A,B>(a, b);
    }
        
    ////////////////////////////////////////////////////////////////////////////////////
    /// Expression template for   A-B
    template<typename A, typename B,
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    class DiffExpression
    {
        const A& a;
        const B& b;

    public:

        typedef typename A::value_type value_type;

        DiffExpression(const A& a, const B& b):a(a), b(b){}
            
        unsigned int size() const { return a.size(); }
            
        const value_type operator[](const unsigned int i) const
        { 
            return a[i] - b[i];
        }
    };
        
    template<typename A, typename B,
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    const DiffExpression<A,B> operator-(const A& a, const B& b)
    {
        return DiffExpression<A,B>(a, b);
    }
        
    ////////////////////////////////////////////////////////////////////////////////////
    /// Expression template for  a*B
    template<typename A, typename B, 
             typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    class LeftScalarMultiplicationExpression
    {
        const A& a;
        const B& b;

    public:

        typedef typename B::value_type value_type;

        LeftScalarMultiplicationExpression(const A& a, const B& b):a(a), b(b){}

        unsigned int size() const { return b.size(); }
            
        const value_type operator[](const unsigned int i) const
        { 
            return a*b[i];
        }
    };
        
    template<typename A, typename B,
             typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    const LeftScalarMultiplicationExpression<A,B> operator*(const A& a, const B& b)
    {
        return LeftScalarMultiplicationExpression<A,B>(a, b);
    }


    ////////////////////////////////////////////////////////////////////////////////////
    /// Expression template for  A*b
    template<typename A, typename B, 
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
    class RightScalarMultiplicationExpression
    {
        const A& a;
        const B& b;

    public:
        typedef typename A::value_type value_type;

        RightScalarMultiplicationExpression(const A& a, const B& b):a(a), b(b){}

        unsigned int size() const { return a.size(); }
            
        const value_type operator[](const unsigned int i) const
        { 
            return a[i]*b;
        }
    };
        
    template<typename A, typename B,
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
    const RightScalarMultiplicationExpression<A,B> operator*(const A& a, const B& b)
    {
        return RightScalarMultiplicationExpression<A,B>(a, b);
    }
    

        
    ////////////////////////////////////////////////////////////////////////////////////
    /// Expression template for  A+b
    template<typename A, typename B,
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
    class RightScalarAdditionExpression
    {
        const A& a;
        const B& b;
    public:
        typedef typename A::value_type value_type;

        RightScalarAdditionExpression(const A& a, const B& b):a(a), b(b){}

        unsigned int size() const { return a.size(); }
            
        const value_type operator[](const unsigned int i) const
        { 
            return a[i]*b;
        }
    };
        
    template<typename A, typename B,
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
    const RightScalarAdditionExpression<A,B> operator+(const A& a, const B& b)
    {
        return RightScalarAdditionExpression<A,B>(a, b);
    }
        
    ////////////////////////////////////////////////////////////////////////////////////
    /// Expression template for  a+B
    template<typename A, typename B,
             typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    class LeftScalarAdditionExpression
    {
        const A& a;
        const B& b;
    public:
        typedef typename B::value_type value_type;

        LeftScalarAdditionExpression(const A& a, const B& b):a(a), b(b){}

        unsigned int size() const { return b.size(); }
            
        const value_type operator[](const unsigned int i) const
        { 
            return b[i]*a;
        }
    };

    template<typename A, typename B,
             typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    const LeftScalarAdditionExpression<A,B> operator+(const A& a, const B& b)
    {
        return LeftScalarAdditionExpression<A,B>(a, b);
    }

    ////////////////////////////////////////////////////////////////////////////////////
    /// Expression template for  A/b
    template<typename A, typename B,
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
    class RightScalarDivisionExpression
    {
        const A& a;
        const B& b;
    public:
        typedef typename A::value_type value_type;

        RightScalarDivisionExpression(const A& a, const B& b):a(a), b(b){}

        unsigned int size() const { return b.size(); }
            
        const value_type operator[](const unsigned int i) const
        { 
            return a[i]/b;
        }
    };

    template<typename A, typename B,
             typename= typename std::enable_if<!std::is_fundamental<A>::value, A>::type,
             typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
    const RightScalarDivisionExpression<A,B> operator/(const A& a, const B& b)
    {
        return RightScalarDivisionExpression<A,B>(a, b);
    }
        
    ////////////////////////////////////////////////////////////////////////////////////
    /// Expression template for  M*A
    template<typename A, typename B, 
             typename= typename std::enable_if<std::is_class<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    class LeftMatrixMultiplicationExpression
    {
        const A& a;
        const B& b;

    public:

        typedef typename A::value_type value_type;

        LeftMatrixMultiplicationExpression(const A& a, const B& b):a(a), b(b){}

        unsigned int size() const { return b.size(); }
            
        const value_type operator[](const unsigned int i) const
        { 
            value_type entry=0;  
            for (index j=0;j<a.shape(0);j++)
                entry+=a.xentry(i,j)*b[j];
            return entry;
        }
    };

    template<typename A, typename B, 
             typename= typename std::enable_if<std::is_class<A>::value, A>::type,
             typename= typename std::enable_if<!std::is_fundamental<B>::value, B>::type>
    const LeftMatrixMultiplicationExpression<A,B> operator*(const A& a, const B& b)
    {
        return LeftMatrixMultiplicationExpression<A,B>(a, b);
    }
}
#endif