//===================================================================
// File:        random.h
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    Pseudo-Random numbers generator.
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


#ifndef _RANDOM_H_
#define _RANDOM_H_


// Supress warnings from external headers
#pragma GCC system_header


#include <ctime>
#include <cassert>
#include <climits>

#include <unistd.h>
#include <sys/time.h>

#include <mpi.h>


#if 0
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#endif


#include <Random123/philox.h>


template <typename R, typename T>
inline R uniform(R val01, T min, T max)
{
	assert(val01 >= 0.0 && val01 < 1.0);

	return val01 * (max - min) + min;
}

template <typename R, typename I>
inline R uniform(I val)
{
	//return val / ((R)std::numeric_limits<I>::max() + 1);
	return val / ((R)UINT_MAX + 1);
}

template <typename R, typename I, typename T>
inline R uniform(I val, T min, T max)
{
	return uniform(uniform<R>(val), min, max);
}


class Seed4
{
private:
	unsigned s[4];
public:
	Seed4()
	{
		int r;
  	timeval tp1;  
  	gettimeofday(&tp1, NULL);

		MPI_Comm_rank(MPI_COMM_WORLD, &r);

		s[0] = r;
		s[1] = getpid();
		s[2] = tp1.tv_sec;
		s[3] = tp1.tv_usec;
	}

	operator const unsigned *() { return s; }
};


#if 0
template <typename Real = double>
class Random
{
	boost::mt19937 rng;
	boost::uniform_01<Real> dist;
public:
	Random(time_t seed = time(NULL)) : rng(seed) {}

	Real operator()() { return dist(rng); }
	Real operator()(Real min, Real max) { return uniform(dist(rng), min, max); }

	operator Real() { return dist(rng); }
	operator unsigned() { return rng(); }
};
#endif


template <typename Real = double>
class Random4x32
{
public:
	enum { N = 4 };
private:
	r123::Philox4x32 rng;
	r123::Philox4x32::key_type key;
	r123::Philox4x32::ctr_type cnt, rnd;
	size_t n;
private:
	void reset()
	{
		cnt.incr();
		rnd = rng(cnt, key);
		n = 0;
	}
	void random()
	{
		if(n >= N) reset();
	}
public:
	Random4x32(unsigned k1, unsigned k2, const unsigned seed[]) : n(N)
	{
		key[0] = k1;
		key[1] = k2;
		
#pragma unroll
		for(int i = 0; i < N; ++i)
			cnt[i] = seed[i];
	}

	Real operator()()
	{
		random();

		assert(n < N);

		return uniform<Real>(rnd[n++]);
	}
	Real operator()(Real min, Real max)
	{
		random();

		assert(n < N);

		return uniform<Real>(rnd[n++], min, max);
	}

	operator Real () { return (*this)(); }
	operator unsigned () { random(); assert(n < N); return rnd[n++]; }
	operator const unsigned * () { reset(); return rnd.v; }
};


#endif /* _RANDOM_H_ */
