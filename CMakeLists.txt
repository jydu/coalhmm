# CMake script for CoalHMM
# Author: Julien Dutheil
# Created: 03/09/2009

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(coalhmm CXX)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF()

SET(CMAKE_CXX_FLAGS "-Wall -Weffc++ -Wshadow")
IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF()

IF(NO_VIRTUAL_COV)
  MESSAGE("-- Covariant return with virtual inheritance disabled.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_VIRTUAL_COV=1")
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

#build info?
IF(NOT DEFINED INFO)
  SET(INFO UNIX)
ENDIF(NOT DEFINED INFO)

#build man pages?
IF(NOT DEFINED MAN)
  SET(MAN UNIX)
ENDIF(NOT DEFINED MAN)

#find executables for documentation
FIND_PROGRAM(MAKEINFO_EXE NAMES makeinfo)
FIND_PROGRAM(PDFTEX_EXE NAMES pdftex)
FIND_PROGRAM(NROFF_EXE NAMES nroff)


#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find build_static)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})

  IF(${build_static})
    SET(${lib_name}_STATIC_NAMES lib${lib_name}.a)
    FIND_LIBRARY(${lib_name}_STATIC_LIBRARY NAMES ${${lib_name}_STATIC_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_STATIC_LIBRARY)
      MESSAGE("-- Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   static libraries: ${${lib_name}_STATIC_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_STATIC_LIBRARY})
  ELSE()
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_LIBRARY)
      MESSAGE("-- Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   dynamic libraries: ${${lib_name}_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
  ENDIF()
  
ENDMACRO(IMPROVED_FIND_LIBRARY)

# Libtool-like version number
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current:  The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age:      The difference between the newest and oldest interfaces that this
#           library implements.
# In other words, the library implements all the interface numbers in the
# range from number current - age to current.
SET(COALHMM_VERSION_CURRENT "1")
SET(COALHMM_VERSION_REVISION "0")
SET(COALHMM_VERSION_AGE "0")

# Effective version number computation
MATH(EXPR COALHMM_VERSION_MAJOR "${COALHMM_VERSION_CURRENT} - ${COALHMM_VERSION_AGE}")
SET(COALHMM_VERSION_MINOR ${COALHMM_VERSION_AGE})
SET(COALHMM_VERSION_PATCH ${COALHMM_VERSION_REVISION})
SET(COALHMM_VERSION "${COALHMM_VERSION_MAJOR}.${COALHMM_VERSION_MINOR}.${COALHMM_VERSION_PATCH}")

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

# Find the dynamic libraries (needed in anycase for building libcoalhmm)
IMPROVED_FIND_LIBRARY(LIBS bpp-phyl Bpp/Phyl/Tree.h 0)
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Bpp/Seq/Alphabet/Alphabet.h 0)
IMPROVED_FIND_LIBRARY(LIBS bpp-core Bpp/Clonable.h 0)

#If a static exe is to be built, then look for static libraries.
#Here the order of check is very important, changing it would lead to a linkage error!
IF(BUILD_STATIC)
  MESSAGE("-- Checking static version of libaries...")
  IMPROVED_FIND_LIBRARY(STATIC_LIBS bpp-phyl Bpp/Phyl/Tree.h 1)
  IMPROVED_FIND_LIBRARY(STATIC_LIBS bpp-seq Bpp/Seq/Alphabet/Alphabet.h 1)
  IMPROVED_FIND_LIBRARY(STATIC_LIBS bpp-core Bpp/Clonable.h 1)
ENDIF(BUILD_STATIC)

IF(MAKEINFO_EXE)
  MESSAGE("-- Found makeinfo here: ${MAKEINFO_EXE}")
  MESSAGE("   Adding targets: info, html")

  IF(INFO)
    SET(ADD_INFO_TO "ALL")
    MESSAGE("   Adding target info to target all")
  ENDIF(INFO)

  ADD_CUSTOM_TARGET(info
    ${ADD_INFO_TO}
    COMMAND ${MAKEINFO_EXE} coalhmm.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )

  ADD_CUSTOM_TARGET(html
    COMMAND ${MAKEINFO_EXE} --html --css-ref=http://www.w3.org/StyleSheets/Core/Steely coalhmm.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
ENDIF(MAKEINFO_EXE)

IF(PDFTEX_EXE)
  MESSAGE("-- Found pdftex here: ${PDFTEX_EXE}")
  MESSAGE("   Adding target: pdf")

  ADD_CUSTOM_TARGET(pdf
    COMMAND ${PDFTEX_EXE} coalhmm.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
ENDIF(PDFTEX_EXE)

IF(NROFF_EXE)
  MESSAGE("-- Found nroff here: ${NROFF_EXE}")
  MESSAGE("   Adding targets: man")

  ADD_CUSTOM_TARGET(man
    ALL
    COMMAND cp coalhmm.1.txt coalhmm.1
    COMMAND gzip -f coalhmm.1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/man
    )

ENDIF(NROFF_EXE)

# Subdirectories
ADD_SUBDIRECTORY(CoalHMM)
ADD_SUBDIRECTORY(Python)
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
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES
 "CMakeFiles"
 "Makefile"
 "_CPack_Packages"
 "CMakeCache.txt"
 ".*\\\\.cmake"
 ".*\\\\.svn"
 ".*\\\\.gz"
 ".*\\\\.zip"
 ".*\\\\.deb"
 ".*\\\\.rpm"
 ".*\\\\..*\\\\.swp"
 ".*stamp"
 "\\\\.sh"
 "CoalHMM/\\\\..*"
 "CoalHMM/coalhmm"
 "doc/coalhmm/"
 "doc/coalhmm\\\\.info"
 "doc/coalhmm\\\\.toc"
 "doc/coalhmm\\\\.vr"
 "doc/coalhmm\\\\.tp"
 "doc/coalhmm\\\\.log"
 "doc/coalhmm\\\\.fn"
 "doc/coalhmm\\\\.ky"
 "doc/coalhmm\\\\.cp"
 "doc/coalhmm\\\\.pg"
 "doc/coalhmm\\\\.aux"
 "doc/coalhmm\\\\.pdf"
 "man/.*\\\\.1.gz"
 "man/.*\\\\.1"
 "debian/coalhmm/"
 "debian/coalhmm\\\\.substvars"
 "debian/coalhmm\\\\.debhelper"
 "debian/debhelper\\\\.log"
 "debian/tmp"
 "\\\\..*\\\\.swp"
 "Examples/\\\\..*"
 "Python/\\\\..*"
 ${CPACK_SOURCE_IGNORE_FILES}
)
IF (MACOS)
  SET(CPACK_GENERATOR "Bundle")
ENDIF()

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
INCLUDE(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
IF (UNIX)
#This creates deb packages:
ADD_CUSTOM_TARGET(deb dpkg-buildpackage)
ADD_DEPENDENCIES(deb dist info man)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info man)

#ADD_DEPENDENCIES(info install)
ENDIF(UNIX)

