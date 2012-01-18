########################################################################
##
## Setup project
##
########################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(iics-v1)
MESSAGE( STATUS "Configuring Interactive Intensive Computing Software" )

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

ADD_EXECUTABLE(v1 common.hpp main.cpp)
TARGET_LINK_YOCTO(v1 cliff mpi)
TARGET_LINK_LIBRARIES(v1 ${MPI_LIBRARIES})
