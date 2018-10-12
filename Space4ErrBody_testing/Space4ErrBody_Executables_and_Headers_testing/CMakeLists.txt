 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

 # Set the header files.
set(SPACE4ERRBODY_HEADERS
  "${CMAKE_CURRENT_SOURCE_DIR}/Space4ErrBody.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/getAlgorithm.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/applicationOutput_tudat.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/applicationOutput_pagmo.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/saveOptimizationResults.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/getStuff.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/updateGuidance.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/updateGuidance_val.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/StopOrNot_val.h"
)

set(SPACE4ERRBODY_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Space4ErrBody.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getCurrentDateTime.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getData.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getFiles.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/updateGuidance.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getNewCoefficients.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/updateGuidance_val.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/StopOrNot_val.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getAngularDistance.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getHeadingToTarget.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getGamma.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getGravs.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getRoot.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/getF.cpp"
)

add_library(Space4ErrBody_Executables_and_Headers_testing STATIC ${SPACE4ERRBODY_SOURCES} ${SPACE4ERRBODY_HEADERS})
setup_tudat_library_target(Space4ErrBody_Executables_and_Headers_testing "${CMAKE_CURRENT_SOURCE_DIR}")
