IF(MPI_FOUND)

ADD_EXECUTABLE(testing-mpi main.cpp
test-bubble.cpp
)

TARGET_LINK_LIBRARIES(testing-mpi sludge)
TARGET_LINK_YOCTO(testing-mpi mpi swamp)
TARGET_LINK_LIBRARIES(testing-mpi ${MPI_LIBRARIES} )

ENDIF(MPI_FOUND)

