###############################################################################
# See 
#   https://cmake.org/Wiki/CMake/Tutorials/How_to_create_a_ProjectConfig.cmake_file
#   https://github.com/forexample/package-example
#   https://cmake.org/Wiki/CMake/Tutorials/Exporting_and_Importing_Targets
#
cmake_minimum_required(VERSION 2.8.12)

PROJECT(NUMCXX C CXX)

set(NUMCXX_MAJOR_VERSION 0)
set(NUMCXX_MINOR_VERSION 11)
set(NUMCXX_PATCH_VERSION 20171113)
set(NUMCXX_VERSION    ${NUMCXX_MAJOR_VERSION}.${NUMCXX_MINOR_VERSION}.${NUMCXX_PATCH_VERSION})

#
# Allow to specify -DNUMCXX_CONFIG=xxx  at cmake time on the 
# command line of cmake.
#
SET(NUMCXX_CONFIG "" CACHE STRING "NUMCXX  Config name")

#
# Define an export name. It seems to be a good idea 
# to that the default is just the package name.
# This export name can be used by find_packge() 
# in the NAMES section in order to select different
# configs registered  for the package
# It should be used for all install(TARGET ... EXPORT ${NUMCXX_EXPORT_NAME}
#
set(NUMCXX_EXPORT_NAME "NUMCXX${NUMCXX_CONFIG}")

#
# Path to a couple of our own FindXXX.cmake files (e.g. UMFPACK)
#
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#################################################################################
#
# Package content
#

set(CMAKE_CXX_STANDARD 11)

find_package(Doxygen)
find_package(LAPACK  REQUIRED QUIET)
find_package(UMFPACK REQUIRED QUIET)

message(STATUS "LAPACK_LIBRARIES : ${LAPACK_LIBRARIES}")
message(STATUS "UMFPACK_LIBRARIES: ${UMFPACK_LIBRARIES}")

#
# LAPACK has no includes... 
#
include_directories(${UMFPACK_INCLUDES})

link_libraries(${LAPACK_LIBRARIES})
link_libraries(${UMFPACK_LIBRARIES})
include_directories(include)


add_subdirectory(src)
add_subdirectory(include)
enable_testing()
add_subdirectory(tests)


#
# Add target dox for building doxygen documentation in  build directory
#
if (DOXYGEN_FOUND)
  configure_file(${PROJECT_SOURCE_DIR}/cmake/Doxyfile.in "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
  # Somehow, doxygen seems to not recognise .md files properly from other
  # than their directory
  add_custom_target(dox 
    COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM )
endif()

#
# Copy numcxx-build to bin. This is a kind of style break
# with regard to the idea of out-of-source build.
#
file(COPY 
  ${PROJECT_SOURCE_DIR}/cmake/numcxx-build
  DESTINATION ${PROJECT_SOURCE_DIR}/bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
  GROUP_READ  GROUP_EXECUTE 
  WORLD_READ WORLD_EXECUTE)


# 
# Archive
# 
set(ZIP_FILE_NAME_ROOT  "numcxx-${NUMCXX_VERSION}")

#
# What is the easy way to do this with CPack proper ?
#
add_custom_target(zip hg archive  ${PROJECT_BINARY_DIR}/${ZIP_FILE_NAME_ROOT}.zip)


#
# Upload doc
#
set(MYWWW "/home/fuhrmann/Wias/www/fuhrmann")
add_custom_target(upload 
                 COMMAND rsync -avu html ${MYWWW}/staging/numcxx
                 COMMAND cp -rp ${PROJECT_BINARY_DIR}/${CPACK_PACKAGE_FILE_NAME}.zip ${MYWWW}/staging/numcxx
                 COMMAND make -C ${MYWWW}/sitesrc upload
                 DEPENDS dox zip
                 )

#######################################################
#######################################################
# From here on, things can be more or less well cloned
# for other packages

########################################################
#
# Magic for package use from the build tree.
#
# This registers the build tree with a global CMake registry,
# On unix, it sits under .cmake/packages. 
# This allows to access the current installation 
# from other packages. No need to "make export" whatsoever.
#
# find_packages()  via the NAMES section is able to distinguish between
# different export names.

export(PACKAGE ${NUMCXX_EXPORT_NAME})

#
# Furthermore, the target export allows client packages to
# figure out if they need a recompile or not.
#
export(TARGETS numcxx 
  FILE "${PROJECT_BINARY_DIR}/${NUMCXX_EXPORT_NAME}Targets.cmake" 
  EXPORT_LINK_INTERFACE_LIBRARIES)


#
# We need to set up  the corresponding 
# Files in ${PROJECT_BINARY_DIR} which is the root of the
# build subdirectory
#
set(NUMCXX_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")
set(NUMCXX_LIBRARIES numcxx)
configure_file(cmake/NUMCXXConfig.cmake.in   "${PROJECT_BINARY_DIR}/${NUMCXX_EXPORT_NAME}Config.cmake" @ONLY)
configure_file(cmake/NUMCXXConfigVersion.cmake.in  "${PROJECT_BINARY_DIR}/${NUMCXX_EXPORT_NAME}ConfigVersion.cmake" @ONLY)



########################################################
#
# Magic for installation 
#
# make install will install the package to ${CMAKE_INSTALL_PREFIX}/{lib|bin|include}
# The cmake find stuff shall be  in lib/cmake/numcxx 
#
# The way it is done here produces a relocatable package.
# https://cmake.org/Wiki/CMake/Tutorials/How_to_create_a_ProjectConfig.cmake_file
#
# It works from the same config.in files 
#  cmake/NUMCXXConfig.cmake.in 
#  cmake/NUMCXXConfigVersion.cmake.in 
# as the export, to be instantiated in another subdirectory in order
# to avoid clashes.


#
# Offer the user the choice of overriding the installation directories
#
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH  "Installation directory for header files")

if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/cmake/numcxx)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH  "Installation directory for CMake files")

# 
# Make relative paths absolute (needed later on)
#
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

#
# Get include dir relative to install directory
#
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"   "${INSTALL_INCLUDE_DIR}")

#
# This will  end up in NUMCXXConfig.cmake making things relocatable
#
set(NUMCXX_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/${REL_INCLUDE_DIR}")

configure_file(cmake/NUMCXXConfig.cmake.in   "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${NUMCXX_EXPORT_NAME}Config.cmake" @ONLY)
configure_file(cmake/NUMCXXConfigVersion.cmake.in  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${NUMCXX_EXPORT_NAME}ConfigVersion.cmake" @ONLY)


#
# Install the ${NUMCXX_EXPORT_NAME}Config.cmake and ${NUMCXX_EXPORT_NAME}ConfigVersion.cmake
#
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${NUMCXX_EXPORT_NAME}Config.cmake"
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${NUMCXX_EXPORT_NAME}ConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" 
  COMPONENT dev)

#
# Install shell script etc.
#
install(PROGRAMS  
  ${PROJECT_SOURCE_DIR}/cmake/numcxx-build
  DESTINATION bin 
  COMPONENT dev)

#
# Install the relocatable export set for use with the install tree
#
# Here, we need to give a relative destination
#
# This creates NUMCXXTargets.cmake and  NUMCXX-<config>.cmake
#
install(EXPORT 
  ${NUMCXX_EXPORT_NAME} 
  DESTINATION  "${DEF_INSTALL_CMAKE_DIR}"
  FILE "${NUMCXX_EXPORT_NAME}Targets.cmake"
  COMPONENT dev)

#
# .. and define an uninstall target
# TODO: rewrite this to be clean CMake
add_custom_target(uninstall cat ${PROJECT_BINARY_DIR}/install_manifest.txt | xargs rm -fv)

