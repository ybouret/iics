
ADD_EXECUTABLE(laponite main.cpp ../config.lua)

TARGET_LINK_LIBRARIES(laponite sludge)
TARGET_LINK_VISIT(laponite)
TARGET_LINK_YOCTO(laponite spade lua++ lua)

INSTALL( TARGETS laponite RUNTIME DESTINATION bin )
