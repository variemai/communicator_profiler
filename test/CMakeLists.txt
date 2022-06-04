cmake_minimum_required(VERSION 3.10)
find_package(MPI REQUIRED)
# set the project name

add_executable(Barrier barrier.c)
target_include_directories(Barrier PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(Barrier PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(Barrier PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(Barrier PUBLIC "$(PROJECT_SOURCE_DIR)/")
target_link_libraries(Barrier PUBLIC mpisee)
