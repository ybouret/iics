########################################################################
##
## Setup project
##
########################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(bubbles++)

########################################################################
##
## using yocto libraries
##
########################################################################
SET( YROOT "$ENV{HOME}/program/yocto4/sdk" )
MESSAGE( STATUS "yocto libraries in ${YROOT}" )
INCLUDE( "${YROOT}/share/yocto/rules.cmake" )

INCLUDE_DIRECTORIES( "${YROOT}/include" )
LINK_DIRECTORIES( "${YROOT}/lib" )

########################################################################
##
## detect hostname and set MPI/VISIT
##
########################################################################
EXECUTE_PROCESS(COMMAND hostname OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
MESSAGE(STATUS "working on <${HOSTNAME}>")

SET(IICS_OK OFF)


## default
IF( NOT IICS_OK )
	FIND_PACKAGE(MPI REQUIRED)
ENDIF()

########################################################################
## main code
########################################################################
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ADD_EXECUTABLE(bubbles common.hpp common.cpp domain.hpp domain.cpp fluid.hpp fluid.cpp main.cpp)
TARGET_LINK_YOCTO(bubbles cliff mpi fs mk)
TARGET_LINK_LIBRARIES(bubbles ${MPI_LIBRARIES})
