CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(laponite-hs)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(SOURCES_SER 
types.hpp 	types.cpp
point.hpp 	point.cpp
spot.hpp  	spot.cpp
bubble.cpp 	bubble.hpp
bubbles.hpp 	bubbles.cpp
)

FIND_PACKAGE(MPI)

IF(MPI_FOUND)
	ADD_DEFINITIONS( -DHAS_MPI )
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	SET(SOURCES_PAR bubble-mpi.cpp)
ENDIF(MPI_FOUND)


ADD_LIBRARY( laphs STATIC
 ${SOURCES_SER}
 ${SOURCES_PAR}
 )

ADD_SUBDIRECTORY(tests)

