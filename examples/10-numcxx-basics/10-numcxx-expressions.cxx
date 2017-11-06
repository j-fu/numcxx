///
/// \example 10-numcxx-expressions.cxx
///
/// Demonstrate the basic usage of arrays and expression templates in numcxx.
///
/// Compile it with
/// ````
/// $ numcxx-compile -o 10-numcxx-expressions 10-numcxx-expressions.cxx
/// ````
/// This invokes cmake to find the installation of the numcxx library,
/// to set up the a small cmake project and then to compile the code.
/// 
///
/// Run with
/// ````
/// $ ./00-hello-world
/// ````
/// 
/// Epression templates  provide two advantages: they allow to work with
/// vectors using the common mathematical notations. Furthermore, for
/// complex expressions, they allow to avoid intermediate storage of results
/// in vectors, providing a possible speed advantage.
///
/// One potential pitfall for expression templates is the use of auto to detect
/// the type of result from an expression template, see e.g. https://eigen.tuxfamily.org/dox/TopicPitfalls.html
/// If you try to do this, please use numcxx::arrayexpr to obtain the proper
/// data type.
#include <iostream>
#include <numcxx/numcxx.hxx>

int main()
{
  const int n=3;
  //////////////////////////////////////////////////////////////
  // Declare three  arrays with double  elements.  Like in  the case
  // with std::vector, the contents of  a class object is placed on
  // the stack, but the large data array is put on the heap.

  numcxx::TArray1<double> A(n);
  numcxx::TArray1<double> B(n);
  numcxx::TArray1<double> C(n);
  numcxx::TArray1<double> D(n);

  //////////////////////////////////////////////////////////////
  // We can assign constant values to an array.
  // The expression template trick allows to do the same for B
  // as the loop does for A.
  for (int i=0;i<n;i++) 
    A(i)=3.0;
  B=3.0;

  std::cout << "Assigning constant values to an array:" << std::endl;
  // numcxx arrays can be printed via iostreams.
  std::cout << "A:"<< std::endl << A << std::endl;
  std::cout << "B:"<< std::endl << B << std::endl;
    
  //////////////////////////////////////////////////////////////
  for (int i=0;i<n;i++) 
    A(i)=i+5;
  B=A;

  std::cout << "Copying data from one array to the other:" << std::endl;
  std::cout << "A:"<< std::endl << A << std::endl;
  std::cout << "B:"<< std::endl << B << std::endl;

  //////////////////////////////////////////////////////////////
  // On arrays of the same size, we can perform the operations
  // allowed for elements of a vector space, i.e. addition, subraction and
  // multiplication by scalars.
  // The expression template trick allows to do the same for C
  // as the loop does for D.
  for (int i=0;i<n;i++)
    C(i)=3*A(i)+B(i);
  D=3*A+B;
  std::cout << "Arithmetic expressions:" << std::endl;
  std::cout << "C:"<< std::endl << C << std::endl;
  std::cout << "D:"<< std::endl << C << std::endl;

   
  //////////////////////////////////////////////////////////////
  // In array expressions, stay safe and avoid 
  // auto on the left hand side!
  numcxx::TArray1<double> E=A+17*B;

  // If you want to use auto, do it like this.
  auto F=numcxx::arrayexpr(A+17*B);

  // The alternative
  // auto F=A+17*B;
  // won't work with << or later use of F on the left
  // hand side of expressions because the type of F detected
  // by auto is the expression template, not the vector.
  // Here, it would  be detected at compile time.
    
  std::cout << "Avoiding auto with expression templates:" << std::endl;
  std::cout << "E:"<< std::endl << E << std::endl;
  std::cout << "F:"<< std::endl << F << std::endl;
    
  //////////////////////////////////////////////////////////////
  // This is the real dangerous situation with auto, 
  // described in https://eigen.tuxfamily.org/dox/TopicPitfalls.html
  auto G=numcxx::arrayexpr(19*A+10*B);
    
  // With
  // auto G=19*A+10*B;
  // G would be an instance of an intermediate type which contains
  // A and B and the possibility to calculate the result of the rhs expression.
  // Changing A would change also A in G and lead to a different result.
  // This is not detected at compile time!
  A(1)++;
  E=G+E;
  std::cout << "Pitfalls of auto:" << std::endl;
  std::cout << "E:"<< std::endl << E << std::endl;
    
}
