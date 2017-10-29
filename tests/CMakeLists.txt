


add_executable(test-expressions test-expressions.cxx)
add_executable(test-jacobi test-jacobi.cxx)
add_executable(test-lapack test-lapack.cxx)
add_executable(test-umfpack test-umfpack.cxx)


add_test(NAME expressions COMMAND  test-expressions WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME jacobi COMMAND  test-jacobi WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME lapack COMMAND  test-lapack WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME umfpack COMMAND  test-umfpack WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


