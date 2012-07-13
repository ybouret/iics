ADD_EXECUTABLE(test-bubbles test-bubbles.cpp)
TARGET_LINK_LIBRARIES(test-bubbles laphs)
TARGET_LINK_YOCTO(test-bubbles)

IF(MPI_FOUND)
	ADD_EXECUTABLE(test-bubbles-mpi test-bubbles-mpi.cpp)
	TARGET_LINK_LIBRARIES(test-bubbles-mpi laphs)
	TARGET_LINK_YOCTO(test-bubbles-mpi mpi)
	TARGET_LINK_LIBRARIES(test-bubbles-mpi ${MPI_LIBRARIES})
ENDIF(MPI_FOUND)


