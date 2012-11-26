CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(sludge2)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../sdk)


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
marker.hpp     marker.cpp
jpack.hpp      jpack.cpp
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
seg/build.cpp
seg/pressurize.cpp
seg/bracket.cpp
seg/io.cpp
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
	mpi/cell.hpp
	mpi/cell.cpp
	mpi/segmenter.cpp
)
SOURCE_GROUP( mpi FILES ${SOURCES_MPI} )

SET(SOURCES_CELL
	mpi/cell/bulk.cpp
	mpi/cell/dispatch.cpp
	mpi/cell/pressure.cpp
	mpi/cell/gradp.cpp
	mpi/cell/velocity.cpp
)
SOURCE_GROUP( mpi/cell FILES ${SOURCES_CELL} )


SET(SOURCES_PAR ${SOURCES_MPI} ${SOURCES_CELL})

ENDIF(MPI_FOUND)

INCLUDE( ${YROOT}/share/yocto/visit.cmake )

IF(VISIT_FOUND)
	SET(SOURCES_VIS_CORE
		visit/simulation.hpp
		visit/simulation.cpp )
	SOURCE_GROUP( visit FILES ${SOURCES_VIS_CORE} )
	
	SET(SOURCES_VIS_SIM
		visit/sim/meta.cpp
		visit/sim/mesh.cpp
		visit/sim/step.cpp
		visit/sim/init.cpp
		visit/sim/data.cpp
		visit/sim/cmds.cpp
		)
	SOURCE_GROUP( visit/sim FILES ${SOURCES_VIS_SIM} )
	
	SET(SOURCES_VIS ${SOURCES_VIS_CORE} ${SOURCES_VIS_SIM})
ENDIF(VISIT_FOUND)

ADD_LIBRARY( sludge STATIC
 ${SOURCES}
 ${SOURCES_PAR}
 ${SOURCES_VIS}
 )

ADD_SUBDIRECTORY(tests)
IF(MPI_FOUND)
	ADD_SUBDIRECTORY(tests-mpi)
ENDIF(MPI_FOUND)

IF(VISIT_FOUND)
	ADD_SUBDIRECTORY(laponite)
ENDIF(VISIT_FOUND)
