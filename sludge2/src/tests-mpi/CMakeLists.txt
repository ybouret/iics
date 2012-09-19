
ADD_EXECUTABLE(testing-mpi 
main.cpp
test-bubbles.cpp
test-param.cpp
)

TARGET_LINK_LIBRARIES(testing-mpi sludge ${MPI_LIBRARIES})
TARGET_LINK_YOCTO(testing-mpi spade mk mpi)

