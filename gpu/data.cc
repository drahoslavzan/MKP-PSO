//===================================================================
// File:        data.cc
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


#include "generic.h"
#include "data.h"

#include <fstream>


MKPData::MKPData(const char *file) throw(const std::string &)
{
	std::ifstream in(file);

	if(!in.good())
		throw std::string("MKPData(): ") + file + ": File doesn't exist";

	in >> n;
	in >> m;
	in >> o;

	if(!in.good())
		throw std::string("MKPData(): ") + file + ": File in wrong format";

	na = (n < 64) ? 64 : nearestPow2(n);

	// Profits
	for(size_t i = 0; i < n; ++i)
	{
		int a;

		if(!in.good())
			throw std::string("MKPData(): ") + file + ": File in wrong format";

		in >> a;

		pro.push_back(a);
	}

	// Align profits
	for(size_t i = 0; i < na - n; ++i)
		pro.push_back(0);

	// Weights
	for(size_t i = 0; i < m; ++i)
	{
		for(size_t j = 0; j < n; ++j)
		{
			int a;

			if(!in.good())
				throw std::string("MKPData(): ") + file + ": File in wrong format";

			in >> a;

			wei.push_back(a);
		}

		// Align weights
		for(size_t j = 0; j < na - n; ++j)
			wei.push_back(0);
	}

	// Capacities
	for(size_t i = 0; i < m; ++i)
	{
		int a;

		if(!in.good())
			throw std::string("MKPData(): ") + file + ": File in wrong format";


		in >> a;

		cap.push_back(a);
	}

	return;
}

