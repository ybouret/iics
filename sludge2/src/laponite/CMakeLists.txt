
ADD_EXECUTABLE(laponite main.cpp)

TARGET_LINK_YOCTO(laponite spade)
TARGET_LINK_LIBRARIES(laponite sludge)
TARGET_LINK_VISIT(laponite)

INSTALL( TARGETS laponite RUNTIME DESTINATION bin )
