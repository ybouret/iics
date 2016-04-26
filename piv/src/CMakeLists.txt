CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(piv)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../)

ADD_EXECUTABLE(genpiv genpiv.cpp sample.lua)

TARGET_LINK_YOCTO(genpiv graphics jpeg png z lua lua++)
INSTALL( TARGETS genpiv RUNTIME DESTINATION bin )
