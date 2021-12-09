cmake_minimum_required(VERSION 3.12)
project(healpeak)

set(CMAKE_CXX_STANDARD 14)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

include_directories(src ${Python3_INCLUDE_DIRS} $ENV{CONDA_PREFIX}/include)
link_directories($ENV{CONDA_PREFIX}/lib)

link_libraries(healpix_cxx cxxsupport sharp cfitsio)

pybind11_add_module(healpeak src/healpeak.cpp)
