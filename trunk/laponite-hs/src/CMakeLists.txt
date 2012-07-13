CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(laponite-hs)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

ADD_LIBRARY( laphs STATIC
types.hpp types.cpp
bubble.cpp bubble.hpp )

ADD_SUBDIRECTORY(tests)

