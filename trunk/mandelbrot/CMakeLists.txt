#introduction
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_SKIP_PREPROCESSED_SOURCE_RULES TRUE) # lighter makefile
SET(CMAKE_SKIP_ASSEMBLY_SOURCE_RULES     TRUE) # lighter makefile
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS    TRUE) # lighter CMakeLists.txt

#a name for the project
PROJECT(mandelbrot)

#for VISIT
EXECUTE_PROCESS(COMMAND hostname OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
MESSAGE(STATUS "working on <${HOSTNAME}>")

SET(VISIT "")


IF( "${HOSTNAME}" MATCHES "math5.*" )
     SET(MPICH2 "/opt/ScaleMP/mpich2/1.4")
     SET(CMAKE_C_COMPILER    "${MPICH2}/bin/mpicc")
     SET(CMAKE_CXX_COMPILER  "${MPICH2}/bin/mpic++")
     MESSAGE( STATUS "flags for math5")
     SET(VISIT "/usr/local/visit2_4_0.linux-x86_64/2.4.0/linux-x86_64")
ELSE()
	SET(VISIT $ENV{VISIT})
ENDIF()

IF( "${VISIT}" STREQUAL "" )
	MESSAGE( FATAL_ERROR "No VisIt" )
ENDIF()


#define default flags
SET(FLAGS "-O2 -Wall")

SET(CMAKE_C_FLAGS ${FLAGS})
SET(CMAKE_CXX_FLAGS ${FLAGS})

INCLUDE_DIRECTORIES(${VISIT}/include/silo/include ${VISIT}/libsim/V2/include)
LINK_DIRECTORIES(${VISIT}/lib ${VISIT}/libsim/V2/lib ${VISIT}/archives)
ADD_EXECUTABLE(mandelbrot mandelbrot.cpp patch.cpp patch.h SimulationExample.h)
TARGET_LINK_LIBRARIES(mandelbrot ${MPI_LIBRARIES} siloh5 simV2 z hdf5)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	TARGET_LINK_LIBRARIES(mandelbrot dl)
ENDIF()

SET(VISIT_UI_DIR "$ENV{HOME}/.visit/ui")
ADD_CUSTOM_COMMAND( TARGET mandelbrot 
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} ARGS -E make_directory  "${VISIT_UI_DIR}"
        COMMAND ${CMAKE_COMMAND} ARGS -E copy
"${CMAKE_CURRENT_SOURCE_DIR}/mandelbrot.ui" "${VISIT_UI_DIR}/"
        COMMENT "Copy mandelbrot.ui to ${VISIT_UI_DIR}"
        )

