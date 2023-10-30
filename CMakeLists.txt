cmake_minimum_required(VERSION 3.20)
project(mpisee VERSION 1.0 DESCRIPTION "mpisee: A communicator-centric profiler")

find_package(MPI REQUIRED CXX)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall")

add_library(mpisee SHARED commprof.cpp utils.cpp)
target_include_directories(mpisee PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(mpisee PRIVATE ${MPI_CXX_LIBRARIES})

if(MPI_CXX_COMPILE_FLAGS)
    target_compile_options(mpisee PRIVATE ${MPI_CXX_COMPILE_FLAGS})
endif()

if(MPI_CXX_LINK_FLAGS)
    set_target_properties(mpisee PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()

set_target_properties(mpisee PROPERTIES SOVERSION 1)
#message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPI_CXX_INCLUDE_DIRS}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

enable_testing()
add_subdirectory(test)
find_program(SRUN_EXEC srun)

if(SRUN_EXEC)
    set(MPIEXEC "${SRUN_EXEC}")
    message(STATUS "Found srun: ${SRUN_EXEC}")
else()
    message(STATUS "srun not found, using default MPIEXEC")
endif()

