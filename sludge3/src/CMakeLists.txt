CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(sludge3)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../sdk)

SET(SOURCES_SEQ
types.hpp   types.cpp
tracer.hpp  tracer.cpp
bubble.hpp  bubble.cpp bubble-ios.cpp bubble-ops.cpp
bubbles.hpp bubbles.cpp
shape.hpp   shape.cpp
)
SOURCE_GROUP( sludge FILES ${SOURCES_SEQ} )

SET(SOURCES_MPI)
FIND_PACKAGE(MPI)
IF(MPI_FOUND)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	MESSAGE( STATUS "include MPI: ${MPI_INCLUDE_PATH}" )
SET(SOURCES_MPI
mpi/bubbles.hpp
mpi/bubbles.cpp
)
SOURCE_GROUP( sludge/mpi FILES ${SOURCES_MPI} )

ENDIF(MPI_FOUND)

ADD_LIBRARY(sludge ${SOURCES_SEQ} ${SOURCES_MPI})
ADD_SUBDIRECTORY(test-seq)

IF(MPI_FOUND)
ADD_SUBDIRECTORY(test-mpi)
ENDIF(MPI_FOUND)
