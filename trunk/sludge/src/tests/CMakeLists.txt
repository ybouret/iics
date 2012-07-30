
ADD_EXECUTABLE(testing main.cpp
test-types.cpp
test-shapes.cpp
test-seg.cpp
test-segpbc.cpp
test-fourier.cpp
)

TARGET_LINK_LIBRARIES(testing sludge)
TARGET_LINK_YOCTO(testing swamp mk)

