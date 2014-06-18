#===================================================================
# File:        CMakeLists.txt
# Author:      Drahoslav Zan
# Email:       izan@fit.vutbr.cz
# Affiliation: Brno University of Technology,
#              Faculty of Information Technology
# Date:        Tue Nov 17 22:12:37 CET 2013
# Comments:    Project build file.
# Project:     A Distributed Multidimensional Knapsack Problem Solver
#              using Island Based Particle Swarm Optimization
#              (DIPSO).
#-------------------------------------------------------------------
# Copyright (C) 2013 Drahoslav Zan
#
# This file is part of DIPSO.
#
# DIPSO is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIPSO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DIPSO. If not, see <http://www.gnu.org/licenses/>.
#===================================================================
# vim: set nowrap sw=2 ts=2


cmake_minimum_required(VERSION 2.6)
project(MKPDIPSO)


set(BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(CMAKE_MODULE_PATH ${BASE_DIR}/libs/cmake)
set(PROJ mkp)


find_package(MPI REQUIRED)
find_package(PAPI REQUIRED)
find_package(CUDA REQUIRED)


set(SOURCES
	data.cc
	main.cc
)

set(CUDASOURCES
	knapsack.cu
)

set(FLAGS
	"-DNDEBUG"
)

set(CFLAGS
	"-W -Wall -Wextra -Wno-long-long -Wno-narrowing"
	"-Wno-unknown-pragmas -Wno-maybe-uninitialized"
	"-Wno-switch -Wno-unused-function"
#	"-pedantic"
	"-O3 -msse4.1 -ffast-math -funroll-loops"
)

set(CUDAFLAGS
	"--pre-include cuda.h -D__STDC_CONSTANT_MACROS"
	"-O3 --use_fast_math"
#	"-m64"
#	"-gencode arch=compute_10,code=sm_10"
	"-gencode arch=compute_20,code=sm_20"
	"-gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35"
)


include_directories(${BASE_DIR}/libs/Random123/include)
include_directories(${MPI_INCLUDE_PATH})

cuda_include_directories(${CUDA_INCLUDE_DIRS})
cuda_include_directories(${CUDA_TOOLKIT_ROOT_DIR}/samples/common/inc)


foreach(arg ${CFLAGS})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${arg}")
endforeach()

foreach(arg ${CUDAFLAGS})
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${arg}")
endforeach()

foreach(arg ${FLAGS})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${arg}")
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${arg}")
endforeach()


cuda_compile(CUDA_O ${CUDASOURCES})

add_executable(${PROJ} ${SOURCES} ${CUDA_O})
target_link_libraries(${PROJ} ${MPI_LIBRARIES} ${CUDA_LIBRARIES} ${PAPI_LIBRARIES})

