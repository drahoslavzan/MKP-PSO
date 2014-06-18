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
#              (MKPDIPSO).
#-------------------------------------------------------------------
# Copyright (C) 2013 Drahoslav Zan
#
# This file is part of MKPDIPSO.
#
# MKPDIPSO is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MKPDIPSO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MKPDIPSO. If not, see <http://www.gnu.org/licenses/>.
#===================================================================
# vim: set nowrap sw=2 ts=2


cmake_minimum_required(VERSION 2.6)
project(MKPDIPSO)


set(BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(CMAKE_MODULE_PATH ${BASE_DIR}/libs/cmake)
set(PROJ mkp)


find_package(MPI REQUIRED)
find_package(PAPI REQUIRED)


set(SOURCES
	main.cc
	knapsack.cc
	data.cc
)

set(CFLAGS
	"-DNDEBUG"
	"-W -Wall -Wextra -Wno-long-long -Wno-narrowing"
	"-Wno-maybe-uninitialized"
	"-Wno-unknown-pragmas"
	"-pedantic -O3"
	"-fopenmp -msse4.1 -ffast-math -funroll-loops"
)


include_directories(${BASE_DIR}/libs/Random123/include)
include_directories(${MPI_INCLUDE_PATH})


foreach(arg ${CFLAGS})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${arg}")
endforeach()


add_executable(${PROJ} ${SOURCES})
target_link_libraries(${PROJ} "-lgomp" ${MPI_LIBRARIES} ${PAPI_LIBRARIES})

