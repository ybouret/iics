CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(fish)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../)

ADD_LIBRARY(fish STATIC fish.hpp fish.cpp
	fish-shell.cpp
	fish-head.cpp
    fish-tail.cpp
	fish-output.cpp
)
ADD_EXECUTABLE(genfish genfish.cpp fish.lua)
TARGET_LINK_LIBRARIES(genfish fish)
TARGET_LINK_YOCTO(genfish mk lua++ lua)

