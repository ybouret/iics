CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(sludge)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(SOURCES_SER 
types.hpp	types.cpp
tracer.hpp	tracer.cpp
bubble.hpp	bubble.cpp
		bubble-vectors.cpp 
		bubble-contour.cpp 
		bubble-shapes.cpp
		bubble-vtk.cpp
spot.hpp	spot.cpp
marker.hpp	marker.cpp
bubbles.hpp	bubbles.cpp
		bubbles-check.cpp
junction.hpp	junction.cpp
segment.hpp	segment.cpp
segmenter.hpp   segmenter.cpp
		segmenter-locate.cpp
		segmenter-assign.cpp
		segmenter-merge.cpp
		segmenter-jpack.cpp
)

FIND_PACKAGE(MPI)

IF(MPI_FOUND)
	ADD_DEFINITIONS( -DHAS_MPI )
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	SET(SOURCES_PAR 
	bubble-mpi.cpp
	bubbles-mpi.cpp
	parameters.hpp	parameters.cpp
	cell.hpp	cell.cpp
			cell-dispatch.cpp
			cell-assemble.cpp
	)
ENDIF(MPI_FOUND)

INCLUDE( ${YROOT}/share/yocto/visit.cmake )

IF( VISIT_FOUND )
	SET(SOURCES_VIS
	)
ENDIF( VISIT_FOUND )

ADD_LIBRARY( sludge STATIC
 ${SOURCES_SER}
 ${SOURCES_PAR}
 ${SOURCES_VIS}
 )

ADD_SUBDIRECTORY(tests)
ADD_SUBDIRECTORY(tests-mpi)


