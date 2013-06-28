
INCLUDE_DIRECTORIES(..)

ADD_EXECUTABLE(laponite main.cpp)

INSTALL( TARGETS laponite RUNTIME DESTINATION bin )
TARGET_LINK_LIBRARIES(laponite sludge)
TARGET_LINK_YOCTO(laponite mk spade lua++ lua)
TARGET_LINK_VISIT(laponite)

