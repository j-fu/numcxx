add_library(numcxx
  simplegrid.cxx
  fem2d.cxx
  fvm2d.cxx
  ../triangle/triangle.c
  )



set_target_properties(numcxx PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(numcxx PROPERTIES COMPILE_FLAGS "-DREAL=double -DTRILIBRARY")
set_property(TARGET numcxx PROPERTY CXX_STANDARD 11)
 
add_executable(triangle  ../triangle/triangle.c)
set_target_properties(triangle PROPERTIES COMPILE_FLAGS "-DREAL=double")
target_link_libraries(triangle m)

add_executable(showme  ../triangle/showme.c)
set_target_properties(showme PROPERTIES COMPILE_FLAGS "-DREAL=double")
target_link_libraries(showme m ${X11_LIBRARIES})

 
install(TARGETS numcxx EXPORT ${NUMCXX_EXPORT_NAME} ARCHIVE DESTINATION ${NUMCXX_REL_INSTALL_DIR}lib)
install(TARGETS triangle showme EXPORT ${NUMCXX_EXPORT_NAME} RUNTIME DESTINATION ${NUMCXX_REL_INSTALL_DIR}bin)


