
INCLUDE_DIRECTORIES(..)

ADD_EXECUTABLE(test-seq main.cpp
test-bubble.cpp
test-grid.cpp
test-curv.cpp
test-segment.cpp
test-ctrl.cpp
)

INSTALL( TARGETS test-seq RUNTIME DESTINATION bin )
TARGET_LINK_LIBRARIES(test-seq sludge)
TARGET_LINK_YOCTO(test-seq mk spade)
