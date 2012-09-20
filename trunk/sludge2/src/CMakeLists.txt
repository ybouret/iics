CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(sludge2)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(SOURCES_SER 
types.hpp      types.cpp
tracer.hpp     tracer.cpp
spot.hpp       spot.cpp
bubble.hpp     bubble.cpp 
bubbles.hpp    bubbles.cpp
junction.hpp   junction.cpp
segment.hpp    segment.cpp
segmenter.hpp  segmenter.cpp 
parameters.hpp parameters.cpp
)

SOURCE_GROUP( core FILES ${SOURCES_SER} )

SET(SOURCES_BUBBLE
bubble/io.cpp
bubble/shapes.cpp
bubble/topology.cpp
)
SOURCE_GROUP( bubble FILES ${SOURCES_BUBBLE} )

SET(SOURCES_SEG
seg/locate.cpp
seg/process.cpp
)
SOURCE_GROUP( seg FILES ${SOURCES_SEG} )


SET(SOURCES ${SOURCES_SER} ${SOURCES_BUBBLE} ${SOURCES_SEG} )

FIND_PACKAGE(MPI)

IF(MPI_FOUND)
MESSAGE( STATUS "Using MPI" )
ADD_DEFINITIONS( -DHAS_MPI )
INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )
SET(SOURCES_MPI
	mpi/bubble.cpp
	mpi/bubbles.cpp
)
SOURCE_GROUP( mpi FILES ${SOURCES_MPI} )
ENDIF(MPI_FOUND)

#INCLUDE( ${YROOT}/share/yocto/visit.cmake )

ADD_LIBRARY( sludge STATIC
 ${SOURCES}
 ${SOURCES_MPI}
 ${SOURCES_VIS}
 )

ADD_SUBDIRECTORY(tests)
IF(MPI_FOUND)
	ADD_SUBDIRECTORY(tests-mpi)
ENDIF(MPI_FOUND)
