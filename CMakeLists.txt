cmake_minimum_required(VERSION 3.4...3.18)
project(compnal)

add_subdirectory(pybind11)
pybind11_add_module(compnal wrapper/src/pybind11_main.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(compnal
#                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})