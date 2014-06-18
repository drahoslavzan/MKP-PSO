//===================================================================
// File:        particle.h
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    The PSO particle.
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


#ifndef _PARTICLE_H_
#define _PARTICLE_H_


#include "bitvector.h"

#include <vector>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <cassert>


class Particle
{
public:
	typedef bitvector<unsigned char> arg_t;
private:
	float best;
	std::vector<float> vel;
	arg_t pos, sol;
public:
	Particle(size_t dimensions)
		: best(-std::numeric_limits<float>::max()),
			vel(dimensions, 0), pos(dimensions, 0), sol(dimensions, 0) {}

	inline float getBest() const
	{
		return best;
	}

	inline void setBest(float v)
	{
		best = v;
		sol = pos;
	}

	inline void setBest(float v, const unsigned char * data)
	{
		best = v;

		//pos.assign(data, data + dims());
		pos = data;
		sol = pos;
	}

	inline float getVelocity(size_t i) const
	{
		assert(i < vel.size());

		return vel[i];
	}

	inline void setVelocity(size_t i, float v)
	{
		assert(i < vel.size());

		vel[i] = v;
	}

	inline float getPosition(size_t i) const
	{
		assert(i < pos.size());

		return pos[i];
	}

	inline void setPosition(size_t i, float v)
	{
		assert(i < pos.size());

		//pos[i] = v;
		pos(i, v);
	}

	inline size_t dims() const
	{
		return pos.size();
	}

	inline const arg_t & solution() const
	{
		return sol;
	}

	inline const arg_t & position() const
	{
		return pos;
	}
};


#endif /* _PARTICLE_H_ */
