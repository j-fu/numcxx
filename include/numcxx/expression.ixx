///
/// \file expression.ixx
///
/// Implementation of expression templates 
///
///
/// Inspired by ideas from:
///
/// Iglberger, K., Hager, G., Treibig,  J., & Rüde, U. (2012). Expression
/// templates    revisited:   a    performance    analysis   of    current
/// methodologies. SIAM Journal on Scientific Computing, 34(2), C42-C69.
/// https://arxiv.org/pdf/1104.1729
///
/// ## Pitfalls encountered.
///
/// ### Pitfall 1
///
/// All classes to be used with the numcxx expression templates
/// should be derived from numcxx specific base classes. The rationale is the control
/// of expression template specialization via ``std::is_base_of``,
/// and the prevention of accidental invocation of the templates
/// in unexpected situations, e.g. with STL vectors.
/// It seems that C++11 provides much better tool to handle this situation via
/// its type_traits
///
/// ### Pitfall 2
///
/// Epression templates using scalars shall store the scalar as a value
/// and not as a reference. E.g. g++ -O optimizes these references away.
/// It seems that this behaviour has to be expected by default,
/// see e.g. Vandevoorde/Josuttis, C++ Templates: The Complete Guide, 2nd Ed. 18.2.1
/// or   http://www.cplusplus.com/forum/general/72582/#msg387184
/// (the broken link therein goes to Vandevoorde/Josuttis) 
///


#ifndef NUMCXX_EXPRESSION_HXX
#define NUMCXX_EXPRESSION_HXX
#include <type_traits>

namespace numcxx
{

  /// Empty base classes for Array expressions
  ///
  /// Used to control template instantiation via type_traits
  class ExpressionBase  { };

  /// Empty base classes for dense matrix aexpressions
  ///
  /// Used to control template instantiation via type_traits
  class MatrixExpressionBase  { };

  /// Empty base classes for sparse matrix expressions
  ///
  /// Used to control template instantiation via type_traits
  class SparseMatrixExpressionBase { };


  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for   A+B
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  class AdditionExpression: public ExpressionBase
  {
    const A& a;
    const B& b;
  public:
    typedef typename A::value_type value_type;
    AdditionExpression(const A& a, const B& b):a(a), b(b){}
    unsigned int size() const { return a.size(); }
    const value_type operator[](const unsigned int i) const {return a[i] + b[i];}
    AdditionExpression& operator=(const AdditionExpression&)=delete;
  };
        
  template<typename A, typename B, 
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  const AdditionExpression<A,B> operator+(const A& a, const B& b)
  {
    return AdditionExpression<A,B>(a, b);
  }


  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for   A-B
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  class SubtractionExpression: public ExpressionBase
  {
    const A& a;
    const B& b;
  public:
    typedef typename A::value_type value_type;
    SubtractionExpression(const A& a, const B& b):a(a), b(b){}
    unsigned int size() const { return a.size(); }
    const value_type operator[](const unsigned int i) const  {   return a[i] - b[i];}
  };
        
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  const SubtractionExpression<A,B> operator-(const A& a, const B& b)
  {
    return SubtractionExpression<A,B>(a, b);
  }
        
  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for  a*B
  template<typename A, typename B, 
           typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  class LeftScalarMultiplicationExpression: public ExpressionBase
  {
    const A a; 
    const B& b;
  public:
    typedef typename B::value_type value_type;
    LeftScalarMultiplicationExpression(const A& a, const B& b):a(a), b(b){}
    unsigned int size() const { return b.size(); }
    const value_type operator[](const unsigned int i) const {return a*b[i];}
  };
        
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  const LeftScalarMultiplicationExpression<A,B> operator*(const A& a, const B& b)
  {
    return LeftScalarMultiplicationExpression<A,B>(a, b);
  }


  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for  A*b
  template<typename A, typename B, 
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
  class RightScalarMultiplicationExpression: public ExpressionBase
  {
    const A& a;
    const B b; 
  public:
    typedef typename A::value_type value_type;
    RightScalarMultiplicationExpression(const A& a, const B& b):a(a),b(b){;}
    unsigned int size() const { return a.size(); }
    const value_type operator[](const unsigned int i) const   {return a[i]*b; }
  };
        
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
  const RightScalarMultiplicationExpression<A,B> operator*(const A& a, const B& b)
  {
    return RightScalarMultiplicationExpression<A,B>(a, b);
  }
    

        
  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for  A+b
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
  class RightScalarAdditionExpression: public ExpressionBase
  {
    const A& a;
    const B b; 
  public:
    typedef typename A::value_type value_type;
    RightScalarAdditionExpression(const A& a, const B& b):a(a), b(b){}
    unsigned int size() const { return a.size(); }
    const value_type operator[](const unsigned int i) const   { return a[i]+b; }
  };
        
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
  const RightScalarAdditionExpression<A,B> operator+(const A& a, const B& b)
  {
    return RightScalarAdditionExpression<A,B>(a, b);
  }
        
  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for  a+B
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  class LeftScalarAdditionExpression: public ExpressionBase
  {
    const A a; 
    const B& b;
  public:
    typedef typename B::value_type value_type;
    LeftScalarAdditionExpression(const A& a, const B& b):a(a), b(b){}
    unsigned int size() const { return b.size(); }
    const value_type operator[](const unsigned int i) const { return a+b[i];}
  };

  template<typename A, typename B,
           typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  const LeftScalarAdditionExpression<A,B> operator+(const A& a, const B& b)
  {
    return LeftScalarAdditionExpression<A,B>(a, b);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for  A-b
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
  class RightScalarSubtractionExpression: public ExpressionBase
  {
    const A& a;
    const B b; 
  public:
    typedef typename A::value_type value_type;
    RightScalarSubtractionExpression(const A& a, const B& b):a(a), b(b){}
    unsigned int size() const { return a.size(); }
    const value_type operator[](const unsigned int i) const   { return a[i]-b; }
  };
        
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
  const RightScalarSubtractionExpression<A,B> operator-(const A& a, const B& b)
  {
    return RightScalarSubtractionExpression<A,B>(a, b);
  }
        
  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for  a-B
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  class LeftScalarSubtractionExpression: public ExpressionBase
  {
    const A a; 
    const B& b;
  public:
    typedef typename B::value_type value_type;
    LeftScalarSubtractionExpression(const A& a, const B& b):a(a), b(b){}
    unsigned int size() const { return b.size(); }
    const value_type operator[](const unsigned int i) const { return a-b[i];}
  };

  template<typename A, typename B,
           typename= typename std::enable_if<std::is_fundamental<A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  const LeftScalarSubtractionExpression<A,B> operator-(const A& a, const B& b)
  {
    return LeftScalarSubtractionExpression<A,B>(a, b);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for  A/b
  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
  class RightScalarDivisionExpression: public ExpressionBase
  {
    const A& a;
    const B b; 
  public:
    typedef typename A::value_type value_type;
    RightScalarDivisionExpression(const A& a, const B& b):a(a), b(b){}
    unsigned int size() const { return a.size(); }
    const value_type operator[](const unsigned int i) const {  return a[i]/b;}
  };

  template<typename A, typename B,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_fundamental<B>::value, B>::type>
  const RightScalarDivisionExpression<A,B> operator/(const A& a, const B& b)
  {
    return RightScalarDivisionExpression<A,B>(a, b);
  }
        
  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for  M*A
  template<typename A, typename B, 
           typename= typename std::enable_if<std::is_base_of<MatrixExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  class LeftMatrixMultiplicationExpression: public ExpressionBase
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
           typename= typename std::enable_if<std::is_base_of<MatrixExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  const LeftMatrixMultiplicationExpression<A,B> operator*(const A& a, const B& b)
  {
    return LeftMatrixMultiplicationExpression<A,B>(a, b);
  }

  ////////////////////////////////////////////////////////////////////////////////////
  /// Expression template for  M*A, M sparse
  template<typename A, typename B, 
           typename= typename std::enable_if<std::is_base_of<SparseMatrixExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  class LeftSparseMatrixMultiplicationExpression: public ExpressionBase
  {
    const A& a;
    const B& b;
  public:
    typedef typename A::value_type value_type;
    LeftSparseMatrixMultiplicationExpression(const A& a, const B& b):a(a), b(b){}
    unsigned int size() const { return b.size(); }
    const value_type operator[](const unsigned int i) const
    { 
      value_type entry=0;  
      auto &IA=*(a.pIA);
      auto &JA=*(a.pJA);
      auto &AA=*(a.pA);

      for (index j=IA[i];j<IA[i+1];j++)
        entry+=AA[j]*b[JA[j]];
      return entry;
    }
  };

  template<typename A, typename B, 
           typename= typename std::enable_if<std::is_base_of<SparseMatrixExpressionBase,A>::value, A>::type,
           typename= typename std::enable_if<std::is_base_of<ExpressionBase,B>::value, B>::type>
  const LeftSparseMatrixMultiplicationExpression<A,B> operator*(const A& a, const B& b)
  {
    return LeftSparseMatrixMultiplicationExpression<A,B>(a, b);
  }
}
#endif
