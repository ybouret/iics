CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(md-bubbles)

SET(YROOT $ENV{HOME}/program/yocto4/sdk)
INCLUDE(${YROOT}/share/yocto/rules.cmake)

INCLUDE_DIRECTORIES( ${YROOT}/include )
LINK_DIRECTORIES(    ${YROOT}/lib     )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../sdk)

ADD_EXECUTABLE(find-bubbles find-bubbles.cpp)
TARGET_LINK_YOCTO(find-bubbles mk)
INSTALL( TARGETS find-bubbles RUNTIME DESTINATION bin )

ADD_EXECUTABLE(parse-data parse-data.cpp)
TARGET_LINK_YOCTO(parse-data)
INSTALL( TARGETS parse-data RUNTIME DESTINATION bin )

ADD_EXECUTABLE(bubble-dyn bubble-dyn.cpp)
TARGET_LINK_YOCTO(bubble-dyn mk)
INSTALL( TARGETS bubble-dyn RUNTIME DESTINATION bin )

