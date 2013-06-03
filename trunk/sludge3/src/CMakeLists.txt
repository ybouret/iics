CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(sludge3)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../sdk)

## Sequential Sources
SET(SOURCES_SEQ
types.hpp     types.cpp
tracer.hpp    tracer.cpp
marker.hpp    marker.cpp
bubble.hpp    bubble.cpp bubble-ios.cpp bubble-ops.cpp bubble-ext.cpp
bubbles.hpp   bubbles.cpp
shape.hpp     shape.cpp
grid.hpp      grid.cpp
junction.hpp  junction.cpp
junctions.hpp junctions.cpp junctions-inter.cpp junctions-segment.cpp
)
SOURCE_GROUP( sludge FILES ${SOURCES_SEQ} )

## MPI Sources
SET(SOURCES_MPI)
FIND_PACKAGE(MPI)
IF(MPI_FOUND)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	MESSAGE( STATUS "include MPI: ${MPI_INCLUDE_PATH}" )
SET(SOURCES_MPI
mpi/bubbles.hpp
mpi/bubbles.cpp
mpi/parameters.hpp
mpi/parameters.cpp
mpi/workspace.hpp
mpi/workspace.cpp
mpi/pressure.cpp mpi/pressure-gradient.cpp
)
SOURCE_GROUP( sludge/mpi FILES ${SOURCES_MPI} )

ENDIF(MPI_FOUND)

##VisIt sources
SET(SOURCES_VIS)
INCLUDE(${YROOT}/share/yocto/visit.cmake)
IF(VISIT_FOUND)
	SET(SOURCES_VIS
	visit/simulation.hpp
	visit/simulation.cpp
	visit/sim-init.cpp
	visit/sim-step.cpp
	)
SOURCE_GROUP( sludge/visit FILES ${SOURCES_VIS} )
ENDIF(VISIT_FOUND)

ADD_LIBRARY(sludge ${SOURCES_SEQ} ${SOURCES_MPI} ${SOURCES_VIS})
ADD_SUBDIRECTORY(test-seq)

IF(MPI_FOUND)
ADD_SUBDIRECTORY(test-mpi)
ENDIF(MPI_FOUND)

IF(VISIT_FOUND)
ADD_SUBDIRECTORY(laponite)
ENDIF(VISIT_FOUND)
