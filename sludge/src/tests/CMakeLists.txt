
ADD_EXECUTABLE(testing main.cpp
test-types.cpp
)

TARGET_LINK_LIBRARIES(testing sludge)
TARGET_LINK_YOCTO(testing)

