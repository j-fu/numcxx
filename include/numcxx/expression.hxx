#ifndef NUMCXX_EXPRESSION_HXX
#define NUMCXX_EXPRESSION_HXX

namespace numcxx
{

    template<typename A, typename B>
    class SumExpression
    {
    public:
        SumExpression(  const A& a, const  B& b ):a_( a ), b_( b ){}
        
        unsigned int size() const { return a_.size(); }
        typedef typename A::value_type value_type;
            
        typename A::value_type xentry( const unsigned int i ) const
        { 
            return a_.xentry(i) + b_.xentry(i);
        }
    private:
        const  A& a_;
        // Reference to the left-hand side operand
        const   B& b_;
        // Reference to the right-hand side operand
    };
        
    template<typename A, typename B>
    const SumExpression<A,B> operator+(const A& a, const B& b )
    {
        return SumExpression<A,B>( a, b );
    }
        
    template<typename A, typename B>
    class DiffExpression
    {
    public:
        DiffExpression(  const A& a, const  B& b ):a_( a ), b_( b ){}
            
        unsigned int size() const { return a_.size(); }
        typedef typename A::value_type value_type;
            
        typename A::value_type xentry( const unsigned int i ) const
        { 
            return a_.xentry(i) - b_.xentry(i);
        }
    private:
        const  A& a_;
        // Reference to the left-hand side operand
        const   B& b_;
        // Reference to the right-hand side operand
    };
        
    template<typename A, typename B>
    const DiffExpression<A,B> operator-(const A& a, const B& b )
    {
        return DiffExpression<A,B>( a, b );
    }
        
        
        
        
    template<typename A, typename B>
    class LeftScalarMultiplicationExpression
    {
    public:
        LeftScalarMultiplicationExpression( const  A& a, const  B& b ):a_( a ), b_( b ){}
        typedef typename B::value_type value_type;
        unsigned int size() const { return b_.size(); }
            
        typename B::value_type xentry( const unsigned int i ) const
        { 
            return a_*b_.xentry(i);
        }
    private:
        const A& a_;
        // Reference to the left-hand side operand
        const  B& b_;
        // Reference to the right-hand side operand
    };
        
    template<typename A, typename B>
    const LeftScalarMultiplicationExpression<A,B> operator*( const A& a, const B& b )
    {
        return LeftScalarMultiplicationExpression<A,B>( a, b );
    }
        
        
        
    template<typename A, typename B>
    class RightScalarDivisionExpression
    {
    public:
        RightScalarDivisionExpression( const  A& a, const  B& b ):a_( a ), b_( b ){}
        typedef typename A::value_type value_type;
        unsigned int size() const { return b_.size(); }
            
        typename A::value_type xentry( const unsigned int i ) const
        { 
            return a_.xentry(i)/b_;
        }
    private:
        const A& a_;
        // Reference to the left-hand side operand
        const  B& b_;
        // Reference to the right-hand side operand
    };
        
    template<typename A, typename B>
    const RightScalarDivisionExpression<A,B> operator/( const A& a, const B& b )
    {
        return RightScalarDivisionExpression<A,B>( a, b );
    }
        
        

}
#endif
