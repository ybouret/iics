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
bubble-io.cpp bubble-contour.cpp bubble-vectors.cpp bubble-shapes.cpp
bubbles.hpp 	bubbles.cpp
intersection.hpp intersection.cpp
segment.hpp	segment.cpp
gmarker.hpp	gmarker.cpp
)

FIND_PACKAGE(MPI)

IF(MPI_FOUND)
	ADD_DEFINITIONS( -DHAS_MPI )
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	SET(SOURCES_PAR 
cell.hpp cell.cpp cell-advect.cpp cell-locate.cpp cell-io.cpp cell-pressure.cpp
bubble-mpi.cpp bubbles-mpi.cpp)
ENDIF(MPI_FOUND)

INCLUDE( ${YROOT}/share/yocto/visit.cmake )

IF( VISIT_FOUND )
	SET(SOURCES_VIS
	simulation.hpp simulation.cpp simulation-visit.cpp)
	
ENDIF( VISIT_FOUND )

ADD_LIBRARY( laphs STATIC
 ${SOURCES_SER}
 ${SOURCES_PAR}
 ${SOURCES_VIS}
 )

ADD_SUBDIRECTORY(tests)

IF( VISIT_FOUND )
	ADD_EXECUTABLE(sim sim.cpp)
	TARGET_LINK_LIBRARIES(sim laphs)
	TARGET_LINK_YOCTO( sim swamp )
	TARGET_LINK_VISIT( sim )
ENDIF( VISIT_FOUND)

