#introduction
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_SKIP_PREPROCESSED_SOURCE_RULES TRUE) # lighter makefile
SET(CMAKE_SKIP_ASSEMBLY_SOURCE_RULES     TRUE) # lighter makefile
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS    TRUE) # lighter CMakeLists.txt

#a name for the project
PROJECT(v0)

#for VISIT
EXECUTE_PROCESS(COMMAND hostname OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
MESSAGE(STATUS "working on <${HOSTNAME}>")

SET(MUST_FIND_MPI ON)

IF( "${HOSTNAME}" MATCHES "math5.*" )
     SET(MUST_FIND_MPI OFF)
     SET(MPICH2 "/opt/ScaleMP/mpich2/1.4")
     SET(CMAKE_C_COMPILER    "${MPICH2}/bin/mpicc")
     SET(CMAKE_CXX_COMPILER  "${MPICH2}/bin/mpic++")
     MESSAGE( STATUS "special flags for math5")
ENDIF()

#VisIt, for everyone, relying on environment
SET(VISIT $ENV{VISIT})
IF( "${VISIT}" STREQUAL "" )
	MESSAGE( FATAL_ERROR "No VisIt" )
ENDIF()

#MPI for everyone but on special architectures
IF(MUST_FIND_MPI)
	FIND_PACKAGE(MPI REQUIRED)
ENDIF()

#define default flags
SET(FLAGS "-O2 -Wall")

SET(CMAKE_C_FLAGS ${FLAGS})
SET(CMAKE_CXX_FLAGS ${FLAGS})

INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${VISIT}/include/silo/include ${VISIT}/libsim/V2/include)
LINK_DIRECTORIES(${VISIT}/lib ${VISIT}/libsim/V2/lib ${VISIT}/archives)
ADD_EXECUTABLE(v0 v0.cpp arrays.c arrays.h types.c types.h v0.ui)
TARGET_LINK_LIBRARIES(v0 ${MPI_LIBRARIES} siloh5 simV2 z hdf5)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	TARGET_LINK_LIBRARIES(v0 dl)
ENDIF()

SET(VISIT_UI_DIR "$ENV{HOME}/.visit/ui")
ADD_CUSTOM_COMMAND( TARGET v0 
	POST_BUILD
	COMMAND ${CMAKE_COMMAND} ARGS -E make_directory  "${VISIT_UI_DIR}"
	COMMAND ${CMAKE_COMMAND} ARGS -E copy "${CMAKE_CURRENT_SOURCE_DIR}/v0.ui" "${VISIT_UI_DIR}/v0.ui"
	COMMENT "Copy v0.ui to ${VISIT_UI_DIR}"
	)
