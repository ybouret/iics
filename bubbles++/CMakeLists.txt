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
ADD_EXECUTABLE(bubbles3d common.hpp common.cpp common3d.hpp common3d.cpp domain3d.hpp domain3d.cpp fluid.hpp fluid.cpp main3d.cpp)
TARGET_LINK_YOCTO(bubbles3d cliff mpi fs mk)
TARGET_LINK_LIBRARIES(bubbles3d ${MPI_LIBRARIES})
