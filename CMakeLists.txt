# CMake script for CoalHMM
# Author: Julien Dutheil
# Created: 03/09/2009

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
PROJECT(coalhmm CXX)

# Compile options
SET (CMAKE_CXX_FLAGS "-std=c++11 -Wall -Weffc++ -Wshadow -Wconversion -Wno-unused-variable -Wno-unused-function")


IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF()


#static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE("-- Static linkage requested. No library will be produced.")
  SET(CMAKE_CXX_FLAGS "-static -static-libgcc ${CMAKE_CXX_FLAGS}")
ENDIF()

# Check compression program
# COMPRESS_PROGRAM controls the choice of program
# COMPRESS_EXT can be used to override the file extension
if (NOT COMPRESS_PROGRAM)
  set (COMPRESS_PROGRAM gzip CACHE STRING "Set program for compressing documentation" FORCE)
endif ()
find_program (COMPRESS_BIN NAMES ${COMPRESS_PROGRAM} DOC "${COMPRESS_PROGRAM} compression program")
if (NOT COMPRESS_BIN)
  message (STATUS "${COMPRESS_PROGRAM} program not found, text doc will not be compressed")
else ()
  # Deduce COMPRESS_EXT for known compression programs if not set
  if (NOT COMPRESS_EXT)
    if (${COMPRESS_PROGRAM} STREQUAL "gzip")
      set (COMPRESS_EXT "gz")
    elseif (${COMPRESS_PROGRAM} STREQUAL "bzip2")
      set (COMPRESS_EXT "bz2")
    else ()
      set (COMPRESS_EXT "${COMPRESS_PROGRAM}") # Default: program name (works for xz/lzma)
    endif ()
  endif ()
  # Generate command line args (always add -c to output compressed file to stdout)
  if (${COMPRESS_PROGRAM} STREQUAL "gzip")
    # -n for no timestamp in files (reproducible builds)
    set (COMPRESS_ARGS -c -n)
  else ()
    set (COMPRESS_ARGS -c)
  endif ()
  message (STATUS "Found ${COMPRESS_BIN} compression program, using file extension .${COMPRESS_EXT}")
endif ()

# Libtool-like version number
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current:  The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age:      The difference between the newest and oldest interfaces that this
#           library implements.
# In other words, the library implements all the interface numbers in the
# range from number current - age to current.
SET(${PROJECT_NAME}_VERSION_CURRENT "1")
SET(${PROJECT_NAME}_VERSION_REVISION "0")
SET(${PROJECT_NAME}_VERSION_AGE "0")

# Effective version number computation
MATH(EXPR ${PROJECT_NAME}_VERSION_MAJOR "${${PROJECT_NAME}_VERSION_CURRENT} - ${${PROJECT_NAME}_VERSION_AGE}")
SET(${PROJECT_NAME}_VERSION_MINOR ${${PROJECT_NAME}_VERSION_AGE})
SET(${PROJECT_NAME}_VERSION_PATCH ${${PROJECT_NAME}_VERSION_REVISION})
SET(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

set (PROJECT_VERSION ${${PROJECT_NAME}_VERSION})

# Find dependencies (add install directory to search)
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

include (GNUInstallDirs)
find_package (bpp-phyl 11.0.0 REQUIRED)

# CMake package
set (cmake-package-location ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
include (CMakePackageConfigHelpers)
configure_package_config_file (
  package.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/package.cmake
  INSTALL_DESTINATION ${cmake-package-location}
  )
write_basic_package_version_file (
  ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  )
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/package.cmake DESTINATION ${cmake-package-location}
  RENAME ${PROJECT_NAME}-config.cmake)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake DESTINATION ${cmake-package-location}
  RENAME ${PROJECT_NAME}-config-version.cmake)

# Subdirectories
add_subdirectory(CoalHMM)
add_subdirectory(Python)
add_subdirectory(doc)
add_subdirectory(man)

# Packager
SET(CPACK_PACKAGE_NAME "coalhmm")
SET(CPACK_PACKAGE_VENDOR "BiRC - Bioinformatics Research Center")
SET(CPACK_PACKAGE_VERSION "1.0.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The CoalHMM program and library")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL")
SET(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
SET(CPACK_SOURCE_IGNORE_FILES 
       "/build/" 
       "/\\\\.git/" 
       "/\\\\.gitignore" 
       ${CPACK_SOURCE_IGNORE_FILES}
       )

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
INCLUDE(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
IF (UNIX)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info man)

#ADD_DEPENDENCIES(info install)
ENDIF(UNIX)

