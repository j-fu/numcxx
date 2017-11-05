cmake_minimum_required(VERSION 2.8.12)

PROJECT(numcxx-examples-numcxx-basics  CXX)

find_package(NUMCXX REQUIRED PATHS /net/wir/lib/CMake/numcxx)
set(CMAKE_CXX_STANDARD 11)

include_directories("${NUMCXX_INCLUDE_DIRS}")

link_libraries("${NUMCXX_LIBRARIES}")

if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")

  add_definitions("-DNUMCXX_CHECK_BOUNDS")

endif()

add_executable(numcxx-expressions  numcxx-expressions.cxx)
add_executable(numcxx-ref  numcxx-ref.cxx)
add_executable(numcxx-sharedptr  numcxx-sharedptr.cxx)

