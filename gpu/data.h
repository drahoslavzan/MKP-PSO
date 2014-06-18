//===================================================================
// File:        data.h
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    Multidimensional Knapsack instance file "parser".
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


#ifndef _DATA_H_
#define _DATA_H_


#include <vector>
#include <cstdlib>
#include <string>


class MKPData
{
	size_t n, m, na;
	int o;
	std::vector<float> pro;
	std::vector<float> wei;
	std::vector<float> cap;
public:
	MKPData(const char *file) throw(const std::string &);

	inline size_t objects() const
	{
		return n;
	}
	
	inline size_t align() const
	{
		return na;
	}
	
	inline size_t knapsacks() const
	{
		return m;
	}
	
	inline int optimum() const
	{
		return o;
	}
	
	inline const float * profits() const
	{
		return pro.data();
	}
	
	inline const float * weights() const
	{
		return wei.data();
	}
	
	inline const float * capacities() const
	{
		return cap.data();
	}
};


#endif /* _DATA_H_ */
