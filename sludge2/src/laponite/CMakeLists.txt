
ADD_EXECUTABLE(laponite main.cpp)

TARGET_LINK_LIBRARIES(laponite sludge)
TARGET_LINK_VISIT(laponite)
TARGET_LINK_YOCTO(laponite spade)

INSTALL( TARGETS laponite RUNTIME DESTINATION bin )
