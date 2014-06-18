//===================================================================
// File:        generic.h
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    Generic components.
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


#ifndef _GENERIC_H_
#define _GENERIC_H_

#include <cstdlib>

#ifdef _OPENMP
#  include <omp.h>
#  define IF_OMP(f) 		f
#  define IF_OMP_0(f) 	f
#  define OMP_PROCS 		omp_get_num_procs()
#else
#  define IF_OMP(f) 		
#  define IF_OMP_0(f) 	0
#  define OMP_PROCS 		1
#endif


#define C_DOUBLE2FLOAT(d) 	d##f
#define C_INT2FLOAT(i) 			i.0f

#define M_DOUBLE2FLOAT(d) 	C_DOUBLE2FLOAT(d)
#define M_INT2FLOAT(i) 			C_INT2FLOAT(i)


#define ELEMS(a) 	(sizeof(a)/sizeof(*a))


#define TO_STR(s) 			#s
#define TO_STR_REF1(r) 	TO_STR(r)


inline size_t nearestPow2(size_t n)
{
	size_t bs = ~0;
	size_t bc = 0;

	for(size_t i = n; i; i >>= 1, bs <<= 1)
		bc += i & 1;

	return (bc > 1) ? (n << 1) & bs : n;
}


// pow(T, unsigned)

template <unsigned n>
struct POW
{
	template <typename T>
	inline static T value(T x) { return x * POW<n - 1>::value(x); }
};

template <>
struct POW<0>
{
	template <typename T>
	inline static T value(T) { return 1; }
};


#endif /* _GENERIC_H_ */
