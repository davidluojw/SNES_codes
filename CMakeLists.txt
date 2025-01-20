cmake_minimum_required(VERSION 3.10)

# Project name and language
project(SNES)

# 指定C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the path to your PETSc installation
set(HOME_DIR /Users/lawkawai) 
set(PETSC_DIR ${HOME_DIR}/lib/petsc-3.22.1-opt)
set(PETSC_ARCH .)

set(CMAKE_C_COMPILER ${HOME_DIR}/lib/mpich-4.2.3-opt/bin/mpicc)
set(CMAKE_CXX_COMPILER ${HOME_DIR}/lib/mpich-4.2.3-opt/bin/mpicxx) 
set(CMAKE_CXX_FLAGS "-O0 -W -Wshadow -Wall -Wextra -Wuninitialized -Wno-unused-parameter")
# set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)

# Find the PETSc package
find_package(PETSc REQUIRED)

include_directories(${PETSC_INC})
include_directories(./include)

add_library(Res_Jac Res_Jac.cpp)
target_link_libraries(Res_Jac ${PETSC_LIB})
# set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
# message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})
# Add an executable for the project
# 添加可执行文件
add_executable(newtons_method SNES_ex1.cpp)
# add_executable(arclength_method SNES_ex16.cpp)
add_executable(test test.cpp)
add_executable(heat_KSP KSP_heat.cpp)
add_executable(heat_SNES SNES_heat.cpp)


# Link PETSc to the project
# 链接PETSc库
target_link_libraries(newtons_method ${PETSC_LIB})
# target_link_libraries(arclength_method ${PETSC_LIB})
target_link_libraries(test ${PETSC_LIB})
target_link_libraries(heat_KSP ${PETSC_LIB})
target_link_libraries(heat_SNES Res_Jac)

# Link Libraries


# Compiler flags and options
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PETSC_C_FLAGS}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PETSC_CXX_FLAGS}")
