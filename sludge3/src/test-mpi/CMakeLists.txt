
INCLUDE_DIRECTORIES(..)

ADD_EXECUTABLE(test-mpi main.cpp
test-bubbles.cpp
test-work.cpp
)

INSTALL( TARGETS test-mpi RUNTIME DESTINATION bin )
TARGET_LINK_LIBRARIES(test-mpi sludge)
TARGET_LINK_YOCTO(test-mpi mk spade mpi)
TARGET_LINK_LIBRARIES(test-mpi ${MPI_LIBRARIES})

