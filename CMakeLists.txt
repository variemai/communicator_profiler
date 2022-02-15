cmake_minimum_required(VERSION 3.10)
project(mpisee VERSION 1.0 DESCRIPTION "MPI Communicator Profiler")

configure_file(commprofConfig.h.in commprofConfig.h)
add_library(mpisee MODULE comprof.c)
