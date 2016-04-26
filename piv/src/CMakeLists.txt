CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(piv)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../)

ADD_EXECUTABLE(genpiv genpiv.cpp)

YOCTO_LINK_TARGET(genpiv graphics jpeg png)
