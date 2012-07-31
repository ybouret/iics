IF(MPI_FOUND)

ADD_EXECUTABLE(testing-mpi main.cpp
test-bubble.cpp
test-bubbles.cpp
test-cell.cpp
test-move.cpp
)

TARGET_LINK_LIBRARIES(testing-mpi sludge)
TARGET_LINK_YOCTO(testing-mpi mpi swamp mk)
TARGET_LINK_LIBRARIES(testing-mpi ${MPI_LIBRARIES} )

ENDIF(MPI_FOUND)

