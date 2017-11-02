add_library(numcxx
  simplegrid.cxx
  ../triangle/triangle.c
  )



set_target_properties(numcxx PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(numcxx PROPERTIES COMPILE_FLAGS "-DREAL=double")
set_property(TARGET numcxx PROPERTY CXX_STANDARD 11)
install(TARGETS numcxx DESTINATION lib EXPORT ${NUMCXX_EXPORT_NAME})


