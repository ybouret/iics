CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(sludge3)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../sdk)

SET(SOURCES_SEQ
types.hpp  types.cpp
tracer.hpp tracer.cpp
)
SOURCE_GROUP( sludge FILES ${SOURCES_SEQ} )

ADD_LIBRARY(sludge ${SOURCES_SEQ})
