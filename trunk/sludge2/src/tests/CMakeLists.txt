
ADD_EXECUTABLE(testing 
main.cpp
test-types.cpp
test-tracers.cpp
test-spots.cpp
test-bubble.cpp
test-seg.cpp
test-param.cpp
test-curv.cpp
)

TARGET_LINK_LIBRARIES(testing sludge)
TARGET_LINK_YOCTO(testing spade mk)

