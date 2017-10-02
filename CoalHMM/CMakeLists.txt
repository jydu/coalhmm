# CMake script for CoalHMM
# Author: Julien Dutheil
# Created: 03/09/2009

SET(coalhmm_H
 CoalHmmStateAlphabet.h
 ThreeSpeciesCoalHmmStateAlphabet.h
 ThreeSpeciesCoalHmmStateAlphabet07.h
 ThreeSpeciesCoalHmmStateAlphabet09.h
 TwoSpeciesDiscretizedCoalHmmStateAlphabet.h
 RateAndCoalHmmStateAlphabet.h
 CoalHmmTransitionMatrix.h
 IlsCoalHmmTransitionMatrix07.h
 IlsCoalHmmTransitionMatrix09.h
 RateAndCoalHmmTransitionMatrix.h
 TwoSpeciesDiscretizedCoalHmmTransitionMatrix.h
 CoalHmmEmissionProbabilities.h
 HomogeneousCoalHmmEmissionProbabilities.h
 HmmTools.h
 CoalHmmSimulator.h
 CtmcTools.h
 generated-two-species/generators.h
 generated-two-species/mapping.h
 NonClockAverageCoalHmmStateAlphabet.h
)

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
  
  # Install bin
  INSTALL(TARGETS coalhmm DESTINATION bin)

ELSE(BUILD_STATIC)
  # Build the static lib
  ADD_LIBRARY(libcoalhmm-static STATIC ${coalhmm_CPP})
  SET_TARGET_PROPERTIES(libcoalhmm-static
    PROPERTIES OUTPUT_NAME coalhmm
    CLEAN_DIRECT_OUTPUT 1
    )
  TARGET_LINK_LIBRARIES(libcoalhmm-static ${LIBS})

  # Build the shared lib
  ADD_LIBRARY(libcoalhmm-shared SHARED ${coalhmm_CPP})
  SET_TARGET_PROPERTIES(libcoalhmm-shared
    PROPERTIES OUTPUT_NAME coalhmm
    CLEAN_DIRECT_OUTPUT 1
    VERSION ${COALHMM_VERSION}
    SOVERSION ${COALHMM_VERSION_MAJOR}
    )
  TARGET_LINK_LIBRARIES(libcoalhmm-shared ${LIBS})
  
  TARGET_LINK_LIBRARIES(coalhmm ${LIBS})

  # Install libs
  INSTALL(TARGETS libcoalhmm-static libcoalhmm-shared DESTINATION lib)

  # Install headers
  INSTALL(FILES ${coalhmm_H} DESTINATION include/CoalHMM)

  # Install bin
  INSTALL(TARGETS coalhmm DESTINATION bin)

ENDIF(BUILD_STATIC)


