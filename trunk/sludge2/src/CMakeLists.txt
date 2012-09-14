CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(sludge2)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(SOURCES_SER 
types.hpp types.cpp
tracer.hpp tracer.cpp
spot.hpp   spot.cpp
)

SOURCE_GROUP( core FILES ${SOURCES_SER} )

INCLUDE( ${YROOT}/share/yocto/visit.cmake )

ADD_LIBRARY( sludge STATIC
 ${SOURCES_SER}
 ${SOURCES_PAR}
 ${SOURCES_VIS}
 )

ADD_SUBDIRECTORY(tests)
