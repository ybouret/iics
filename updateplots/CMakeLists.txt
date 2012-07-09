#introduction
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_SKIP_PREPROCESSED_SOURCE_RULES TRUE) # lighter makefile
SET(CMAKE_SKIP_ASSEMBLY_SOURCE_RULES     TRUE) # lighter makefile
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS    TRUE) # lighter CMakeLists.txt

#a name for the project
PROJECT(updateplots)

#for VISIT
EXECUTE_PROCESS(COMMAND hostname OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
MESSAGE(STATUS "working on <${HOSTNAME}>")

SET(VISIT "")

IF( "${HOSTNAME}" MATCHES "math7.*" )
# do not forget to source /opt/intel/Compiler/11.1/072/bin/iccvars.sh intel64
	FIND_PACKAGE(MPI REQUIRED)
	MESSAGE( STATUS "flags for math7")
	SET(VISIT "/usr/local/visit/2.4.0/linux-x86_64")
ENDIF()

#IF( "${HOSTNAME}" MATCHES "math5.*" )
#	#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/OpenMPI/lib
#	#and launch with '/usr/local/OpenMPI/bin/mpirun -np ...'
#	SET(OPENMPI "/usr/local/OpenMPI")
#	SET(CMAKE_C_COMPILER    "${OPENMPI}/bin/mpicc")
#	SET(CMAKE_CXX_COMPILER  "${OPENMPI}/bin/mpic++")
#       MESSAGE( STATUS "flags for math5")
#        SET(VISIT "/usr/local/visit2_4_0.linux-x86_64/2.4.0/linux-x86_64")
#ENDIF()


IF( "${HOSTNAME}" MATCHES "math5.*" )
     SET(MPICH2 "/opt/ScaleMP/mpich2/1.4")
     SET(CMAKE_C_COMPILER    "${MPICH2}/bin/mpicc")
     SET(CMAKE_CXX_COMPILER  "${MPICH2}/bin/mpic++")
     MESSAGE( STATUS "flags for math5")
     SET(VISIT "/usr/local/visit2_4_0.linux-x86_64/2.4.0/linux-x86_64")
ENDIF()

IF( "${HOSTNAME}" MATCHES "pcmath162.*" )
	FIND_PACKAGE(MPI REQUIRED)
        MESSAGE( STATUS "flags for pcmath162")
        SET(VISIT "/usr/local/visit/2.4.0/linux-x86_64")
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	FIND_PACKAGE(MPI REQUIRED)
	MESSAGE( STATUS "flags for Mac Os X")
#	SET(VISIT "/usr/local/visit/2.4.0/darwin-x86_64") 
	SET(VISIT "/Applications/VisIt.app/Contents/Resources/2.4.1/darwin-x86_64") 
ENDIF()

#define default flags
SET(FLAGS "-O2 -g")

SET(CMAKE_C_FLAGS ${FLAGS})
SET(CMAKE_CXX_FLAGS ${FLAGS})

INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${VISIT}/include/silo/include ${VISIT}/libsim/V2/include)
LINK_DIRECTORIES(${VISIT}/lib ${VISIT}/libsim/V2/lib ${VISIT}/archives)
ADD_EXECUTABLE(updateplots updateplots.cpp SimulationExample.h updateplots.ui)
TARGET_LINK_LIBRARIES(updateplots ${MPI_LIBRARIES} siloh5 simV2 z hdf5)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	TARGET_LINK_LIBRARIES(updateplots dl)
ENDIF()

SET(VISIT_UI_DIR "$ENV{HOME}/.visit/ui")
ADD_CUSTOM_COMMAND( TARGET updateplots 
	POST_BUILD
	COMMAND ${CMAKE_COMMAND} ARGS -E make_directory  "${VISIT_UI_DIR}"
	COMMAND ${CMAKE_COMMAND} ARGS -E copy "${CMAKE_CURRENT_SOURCE_DIR}/updateplots.ui" "${VISIT_UI_DIR}/updateplots.ui"
	COMMENT "Copy updateplots.ui to ${VISIT_UI_DIR}"
	)
