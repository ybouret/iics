########################################################################
##
## Setup project
##
########################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(vdw)

########################################################################
##
## using yocto libraries
##
########################################################################
SET( YROOT "$ENV{HOME}/program/yocto4/sdk" )
MESSAGE( STATUS "yocto libraries in ${YROOT}" )
INCLUDE( "${YROOT}/share/yocto/rules.cmake" )

INCLUDE_DIRECTORIES( "${YROOT}/include" )
LINK_DIRECTORIES( "${YROOT}/lib" )



ADD_EXECUTABLE(vdw main.cpp)
TARGET_LINK_YOCTO(vdw mk)
