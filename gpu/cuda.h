//===================================================================
// File:        cuda.h
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    The CUDA headers.
// Project:     A Distributed Multidimensional Knapsack Problem Solver
//              using Island Based Particle Swarm Optimization
//              (MKPDIPSO).
//-------------------------------------------------------------------
// Copyright (C) 2013 Drahoslav Zan
//
// This file is part of MKPDIPSO.
//
// MKPDIPSO is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MKPDIPSO is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MKPDIPSO. If not, see <http://www.gnu.org/licenses/>.
//===================================================================
// vim: set nowrap sw=2 ts=2


#ifndef _CUDA_H_
#define _CUDA_H_


#ifdef __CUDACC__
#  undef _GLIBCXX_ATOMIC_BUILTINS
#  undef _GLIBCXX_USE_INT128
#endif


#ifdef __CUDACC__
#  include <cassert>
#  include <cfloat>
#  include <cuda_runtime.h>
#  include <cuda_runtime_api.h>
#  include <helper_cuda.h>
#  include <helper_functions.h>
#endif


#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#  undef  assert
#  define assert(arg)
#endif


#define KERNEL 						__global__ void
#define DEVICE(kern, n) 	kern<<<n / 128 + 1, 128>>>


#endif /* _CUDA_H_ */
