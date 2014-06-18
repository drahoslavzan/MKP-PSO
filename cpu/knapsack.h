//===================================================================
// File:        knapsack.h
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    MKP problem solver.
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


#ifndef _KNAPSACK_H_
#define _KNAPSACK_H_


#include "particle.h"
#include "data.h"

#include <vector>
#include <cstring>

#include <mpi.h>


// ===================================================================
// PSO parameters from
//   A Novel Probability Binary Particle Swarm Optimization Algorithm
//   and Its Application.
//   JOURNAL OF SOFTWARE, VOL. 3, NO. 9, DECEMBER 2008
// ===================================================================

const size_t PSO_MIN_PARTICLES = 4;

const size_t PSO_DEFAULT_EXCHANGE_INTERVAL        = 1;
const size_t PSO_DEFAULT_STATS_INTERVAL           = PSO_DEFAULT_EXCHANGE_INTERVAL;

const size_t PSO_DEFAULT_EPOCHS                   = 2000;
const size_t PSO_DEFAULT_STOP_CRITERION           = 0;
const size_t PSO_DEFAULT_PARTICLES                = 1 * PSO_MIN_PARTICLES;

const float PSO_DEFAULT_INERTIA_WEIGHT            = 1.f;
const float PSO_DEFAULT_COGNITIVE                 = 2.f;
const float PSO_DEFAULT_SOCIAL                    = PSO_DEFAULT_COGNITIVE;

const float MKP_DEFAULT_PENALTY                   = 5000.f;
const float MKP_DEFAULT_VELOCITY_CLAMP            = 10.f;

const int MPIROOT = 0;


class MKPSolver
{
public:
	typedef Particle::arg_t sol_t;
public:
	struct Stats
	{
		struct
		{
			float best, worst, avg;
		} fitness;

		std::vector<size_t> weight;
		size_t iterations;
	};
private:
	const MKPData & data;
	std::vector<Particle> swarm;
	float penalty;
	size_t dims;
	size_t best;
	float vMax;
	float iW, lC1, lC2;
private:
	int procs, id;
	float * mstats;
private:
	float exLastFitness;
	Particle::arg_t::elem_t * hsols;
private:
	mutable Stats stats;
private:
	void init();
	bool findBest();
	void update();
	void fitness(const Particle::arg_t * x[4], float * v);
private:
#ifdef SHOW_SYNC
	size_t showStats(size_t n, bool sync = true);
#else
	size_t showStats(size_t n);
#endif
	bool exchangeSync(size_t worst);
public:
	static size_t adjustParticleCount(size_t n);
	static size_t adjustEpochs(size_t n);
public:
	MKPSolver(const MKPData &dat, size_t particles = PSO_DEFAULT_PARTICLES);
	~MKPSolver();

	void setInertiaWeight(float w = PSO_DEFAULT_INERTIA_WEIGHT);
	void setLearningFactors(float c1 = PSO_DEFAULT_COGNITIVE,
			float c2 = PSO_DEFAULT_SOCIAL);
	void setVelocityMax(float m = MKP_DEFAULT_VELOCITY_CLAMP);
	void setPenalty(float p = MKP_DEFAULT_PENALTY);

	const sol_t & solve(size_t epochs = PSO_DEFAULT_EPOCHS,
			size_t stop = PSO_DEFAULT_STOP_CRITERION,
			size_t exchange = PSO_DEFAULT_EXCHANGE_INTERVAL,
			size_t statistics = PSO_DEFAULT_STATS_INTERVAL);

	const Stats & statistics() const;
};


#endif /* _KNAPSACK_H_ */
