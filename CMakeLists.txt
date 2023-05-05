cmake_minimum_required(VERSION 3.15)
project(healpeak)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

include_directories(src ${Python3_INCLUDE_DIRS} $ENV{CONDA_PREFIX}/include)
link_directories($ENV{CONDA_PREFIX}/lib)

link_libraries(healpix_cxx cxxsupport sharp cfitsio)

pybind11_add_module(healpeak src/healpeak.cpp)
