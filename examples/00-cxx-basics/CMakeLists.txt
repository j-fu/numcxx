cmake_minimum_required(VERSION 2.8.12)
project(numcxx-examples-cxx-basics C CXX)
set(CMAKE_CXX_STANDARD 11)

add_executable(01-c-style-stack EXCLUDE_FROM_ALL 01-c-style-stack.cxx)
add_executable(02-c-style-heap  EXCLUDE_FROM_ALL 02-c-style-heap.cxx)
add_executable(03-cxx-style-ref EXCLUDE_FROM_ALL 03-cxx-style-ref.cxx)
add_executable(04-cxx-style-sharedptr EXCLUDE_FROM_ALL 04-cxx-style-sharedptr.cxx)
add_executable(05-cxx-style-sharedptr2 EXCLUDE_FROM_ALL 05-cxx-style-sharedptr2.cxx)



