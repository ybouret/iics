
INCLUDE_DIRECTORIES(..)

ADD_EXECUTABLE(test-seq main.cpp
test-bubble.cpp)

INSTALL( TARGETS test-seq RUNTIME DESTINATION bin )
TARGET_LINK_LIBRARIES(test-seq sludge)
TARGET_LINK_YOCTO(test-seq mk)
