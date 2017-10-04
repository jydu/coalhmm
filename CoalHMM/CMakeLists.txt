# CMake script for CoalHMM
# Author: Julien Dutheil
# Created: 03/09/2009

SET(coalhmm_CPP
 ThreeSpeciesCoalHmmStateAlphabet07.cpp
 ThreeSpeciesCoalHmmStateAlphabet09.cpp
 TwoSpeciesWithOutgroupDiscretizedCoalHmmStateAlphabet.cpp
 TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet.cpp
 RateAndCoalHmmStateAlphabet.cpp
 CoalHmmEmissionProbabilities.cpp
 HomogeneousCoalHmmEmissionProbabilities.cpp
 IlsCoalHmmTransitionMatrix07.cpp
 IlsCoalHmmTransitionMatrix09.cpp
 RateAndCoalHmmTransitionMatrix.cpp
 TwoSpeciesDiscretizedCoalHmmTransitionMatrix.cpp
 HmmTools.cpp
 CoalHmmSimulator.cpp
 CtmcTools.cpp
 generated-two-species/generators.cpp
 generated-two-species/mapping.cpp
 generated-two-species/statesets.cpp
 NonClockAverageCoalHmmStateAlphabet.cpp
)

# Build program
ADD_EXECUTABLE(coalhmm ${coalhmm_CPP} CoalHMM.cpp)

IF(BUILD_STATIC)
  TARGET_LINK_LIBRARIES(coalhmm ${STATIC_LIBS})
  SET_TARGET_PROPERTIES (coalhmm PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  
ELSE(BUILD_STATIC)
  # Build the static lib
  ADD_LIBRARY(libcoalhmm-static STATIC ${coalhmm_CPP})
  target_include_directories (libcoalhmm-static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  )
  set_target_properties (libcoalhmm-static PROPERTIES OUTPUT_NAME libcoalhmm)
  target_link_libraries (libcoalhmm-static ${BPP_LIBS_STATIC})


  # Build the shared lib
  ADD_LIBRARY(libcoalhmm-shared SHARED ${coalhmm_CPP})
  target_include_directories (libcoalhmm-shared PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  )
  set_target_properties (libcoalhmm-shared
    PROPERTIES OUTPUT_NAME libcoalhmm
    MACOSX_RPATH 1
    VERSION ${${PROJECT_NAME}_VERSION}
    SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}
  )
  target_link_libraries (libcoalhmm-shared ${BPP_LIBS_SHARED})


  # Link
  TARGET_LINK_LIBRARIES (coalhmm ${BPP_LIBS_SHARED})
  SET_TARGET_PROPERTIES (coalhmm PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

  # Install libs and headers
  install (
    TARGETS libcoalhmm-static libcoalhmm-shared
    EXPORT libcoalhmm-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
  install (
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/CoalHMM
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
  )
  # Install cmake file representing targets
  install (EXPORT libcoalhmm-targets DESTINATION ${cmake-package-location})
  
ENDIF(BUILD_STATIC)
  
# Install bin
INSTALL(TARGETS coalhmm DESTINATION ${CMAKE_INSTALL_BINDIR})



