ADD_EXECUTABLE(test-bubble test-bubble.cpp)
TARGET_LINK_LIBRARIES(test-bubble laphs)
TARGET_LINK_YOCTO(test-bubble)

ADD_EXECUTABLE(test-pbc test-pbc.cpp)
TARGET_LINK_LIBRARIES(test-pbc laphs yocto)


IF(MPI_FOUND)
	ADD_EXECUTABLE(test-bubble-mpi test-bubble-mpi.cpp)
        TARGET_LINK_LIBRARIES(test-bubble-mpi laphs)
        TARGET_LINK_YOCTO(test-bubble-mpi mpi)
        TARGET_LINK_LIBRARIES(test-bubble-mpi ${MPI_LIBRARIES})

	ADD_EXECUTABLE(test-bubbles-mpi test-bubbles-mpi.cpp)
	TARGET_LINK_LIBRARIES(test-bubbles-mpi laphs)
	TARGET_LINK_YOCTO(test-bubbles-mpi mpi)
	TARGET_LINK_LIBRARIES(test-bubbles-mpi ${MPI_LIBRARIES})

	ADD_EXECUTABLE(test-cell test-cell.cpp)
	TARGET_LINK_LIBRARIES(test-cell laphs)
	TARGET_LINK_YOCTO(test-cell mpi swamp)
        TARGET_LINK_LIBRARIES(test-cell ${MPI_LIBRARIES})
ENDIF(MPI_FOUND)


