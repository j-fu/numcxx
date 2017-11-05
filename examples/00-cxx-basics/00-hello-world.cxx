/**
   \example 00-hello-world.cxx
   
   My first program in C++.

   Compile it with
   ````
   $ g++ -o 00-hello-world 00-hello-world.cxx
   ````
   This calls the compiler ``g++``, tells it to
   compile ``00-hello-world.cxx`` and to create an executable
   named ``00-hello-world``

   Run with
   ````
   $ ./00-hello-world
   ````

 */




// In order to use the input/output library, we need to 
// include a header which describes its interface.
#include <iostream>


// Every C++ program starts with a function called main()
// Invoking the executable 
int main()
{
  // std::cout is the output stream of the program.
  // The << operator writes data to it.
  // std::endl outputs an end-of-line symbol.
  std::cout << "Hello World" << std::endl;
}
