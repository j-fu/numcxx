add_library(numcxx
  simplegrid.cxx
  ../triangle/triangle.c
  )



set_target_properties(numcxx PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(numcxx PROPERTIES COMPILE_FLAGS "-DREAL=double -DTRILIBRARY")
set_property(TARGET numcxx PROPERTY CXX_STANDARD 11)
  
install(TARGETS numcxx EXPORT ${NUMCXX_EXPORT_NAME} ARCHIVE DESTINATION ${NUMCXX_REL_INSTALL_DIR}lib)


