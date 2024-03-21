cmake_minimum_required(VERSION 3.10)
project(mpisee VERSION 1.0 DESCRIPTION "mpisee: A communicator-centric profiler")

find_package(MPI REQUIRED CXX)
find_package(SQLite3 REQUIRED)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall")
# Set default values
set(NUM_BUCKETS 8 CACHE STRING "Number of buckets")
set(BUCKETS "128,1024,8192,65536,262144,1048576,33554432" CACHE STRING "Message size buckets in bytes")

# Generate the header file (the configure_file part remains the same)
configure_file(utils.h.in ${CMAKE_SOURCE_DIR}/utils.h)


add_library(mpisee SHARED commprof.cpp utils.cpp create_db.cpp)
#add_executable(create_db create_db.cpp utils.cpp)
#target_include_directories(create_db PRIVATE ${MPI_CXX_INCLUDE_DIRS})
#target_link_libraries(create_db PRIVATE ${MPI_CXX_LIBRARIES})
target_include_directories(mpisee PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(mpisee PRIVATE ${MPI_CXX_LIBRARIES})
target_link_libraries(mpisee PRIVATE ${SQLite3_LIBRARIES})
target_include_directories(mpisee PRIVATE ${SQLite3_INCLUDE_DIRS})

if(MPI_CXX_COMPILE_FLAGS)
    target_compile_options(mpisee PRIVATE ${MPI_CXX_COMPILE_FLAGS})
endif()

if(MPI_CXX_LINK_FLAGS)
    set_target_properties(mpisee PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()

set_target_properties(mpisee PROPERTIES SOVERSION 1)
#message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPI_CXX_INCLUDE_DIRS}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

#enable_testing()
#add_subdirectory(test)
find_program(SRUN_EXEC srun)

if(SRUN_EXEC)
    set(MPIEXEC "${SRUN_EXEC}")
    message(STATUS "Found srun: ${SRUN_EXEC}")
else()
    message(STATUS "srun not found, using default MPIEXEC")
endif()

