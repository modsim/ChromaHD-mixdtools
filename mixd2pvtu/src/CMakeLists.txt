cmake_minimum_required(VERSION 2.8)
PROJECT(mixd2pvtu)

# Require MPI for this project:

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
# add_compile_options(-g -ftrapv)
add_compile_options(-O3 -Ofast -march=native -mtune=native)

find_package(VTK REQUIRED)
# include(${VTK_USE_FILE})

# set(VTK_DIR "/usr/lib64/vtk")
# set(LDFLAGS "-Wl,--copy-dt-needed-entries")
# set(VTK_DIR "/usr/include/vtk-7.1")

set(CMAKE_INSTALL_PREFIX "~/local")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
link_libraries(stdc++fs)

include_directories(.)
file(GLOB SOURCES "*.cpp")

add_executable(mixd2pvtu ${SOURCES})

target_link_libraries(mixd2pvtu ${MPI_C_LIBRARIES})
target_link_libraries(mixd2pvtu ${MPI_CXX_LIBRARIES})
target_link_libraries(mixd2pvtu ${VTK_LIBRARIES})


install (TARGETS mixd2pvtu DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
