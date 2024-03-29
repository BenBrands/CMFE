CMAKE_MINIMUM_REQUIRED( VERSION 3.1 FATAL_ERROR )

FIND_PACKAGE( deal.II 9.0.0 REQUIRED 
    HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
    )
IF( NOT ${deal.II_FOUND} )
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

IF(NOT (DEAL_II_WITH_PETSC OR DEAL_II_WITH_TRILINOS) OR NOT DEAL_II_WITH_P4EST)
    MESSAGE(FATAL_ERROR "
    Error! The deal.II library found at ${DEAL_II_PATH} was not configured with 
        DEAL_II_WITH_PETSC = ON 
        DEAL_II_WITH_P4EST = ON 
    or 
        DEAL_II_WITH_TRILINOS = ON 
        DEAL_II_WITH_P4EST = 0N
    One or both of these combinations are OFF in your installation but at least one is required for this project."
        )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

INCLUDE( ${DEAL_II_FEATURE_CONFIG} )

PROJECT( ContinuumMechanics )

## find the Boost libraries
set( BOOST_ROOT ${DEAL_II_BOOST_DIR} )
find_package( Boost REQUIRED COMPONENTS filesystem serialization mpi )
INCLUDE_DIRECTORIES( SYSTEM ${Boost_INCLUDE_DIR} )
LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS} )

## Binaries will be saved in given directory
SET(PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}/bin)
## executables will be saved in given directory
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_subdirectory( Executables )
add_subdirectory( ConstitutiveLaws )
add_subdirectory( FE )

ENABLE_TESTING()
ADD_SUBDIRECTORY(tests/ConvergenceTests)