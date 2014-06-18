//===================================================================
// File:        knapsack.cc
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Thu Dec 19 15:51:17 CET 2013
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


#include "knapsack.h"
#include "generic.h"
#include "random.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#include <mpi.h>


using namespace std;


inline static float clamp(float v, float m)
{
	if(v < -m) return -m;
	if(v > m) return m;
	return v;
}

inline static float poslin(float v)
{
	return (v > 0) ? v : 0;
}


size_t MKPSolver::adjustParticleCount(size_t n)
{
	const size_t N = PSO_MIN_PARTICLES;

	return (n % N) ? n + (N - (n % N)) : n;
}

size_t MKPSolver::adjustEpochs(size_t n)
{
	return n;
}


void MKPSolver::init()
{
#pragma omp parallel
	{
		Random4x32<float> rnd(IF_OMP_0(omp_get_thread_num()), 0xcaffe123, Seed4());

#pragma omp for schedule(static)
		for(size_t i = 0; i < swarm.size(); ++i)
		{
			for(size_t j = 0; j < dims; ++j)
			{
				swarm[i].setPosition(j, rnd() < 0.5f ? 1 : 0);
				swarm[i].setVelocity(j, rnd() * vMax);
			}
		}
	}

	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	hsols = new Particle::arg_t::elem_t[dims];

	size_t mstatssz;

	if(id == MPIROOT) mstatssz = 3 * procs;
	else mstatssz = 3;

	mstats = new float[mstatssz];
	memset(mstats, 0, mstatssz);
}

bool MKPSolver::findBest()
{
	bool nb = false;

	// 4x unroll
#pragma omp for schedule(static)
	for(size_t i = 0; i < swarm.size(); i += 4)
	{
		float fval[4];
		const Particle::arg_t * x[] = {
			&swarm[i + 0].position(), &swarm[i + 1].position(),
			&swarm[i + 2].position(), &swarm[i + 3].position() };

		fitness(x, fval);

		float fbest = swarm[best].getBest();
		size_t ibest;

#pragma unroll
		for(int j = 0; j < 4; ++j)
		{
			if(fval[j] > swarm[i + j].getBest())
			{
				swarm[i + j].setBest(fval[j]);

				if(fval[j] > fbest)
				{
					fbest = fval[j];
					ibest = i + j;
				}
			}
		}

		if(fbest > swarm[best].getBest())
		{
#pragma omp critical
			if(fbest > swarm[best].getBest())
			{
				best = ibest;
				nb = true;
			}
		}
	}

	return nb;
}

void MKPSolver::update()
{
	Random4x32<float> rnds(IF_OMP_0(omp_get_thread_num()), 0xabcf0198, Seed4());

	// 4x unroll

	assert(Random4x32<>::N == 4);
	assert(swarm.size() % 4 == 0);

#pragma omp for schedule(static)
	for(size_t i = 0; i < swarm.size(); i += 4)
	{
		const Particle::arg_t * pbest[] = {
			&swarm[i + 0].solution(), &swarm[i + 1].solution(),
			&swarm[i + 2].solution(), &swarm[i + 3].solution()
		};

		const Particle::arg_t & gbest = swarm[best].solution();

		for(size_t n = 0; n < 4; ++n)
			assert(pbest[n]->size() == dims);

		assert(gbest.size() == dims);

		const Particle::arg_t * pos[] = {
			&swarm[i + 0].position(), &swarm[i + 1].position(),
			&swarm[i + 2].position(), &swarm[i + 3].position()
		};

		for(size_t j = 0; j < dims; ++j)
		{
			const unsigned *rnd = rnds;

			float v[4];

#pragma unroll
			for(size_t n = 0; n < 4; ++n)
				v[n] = iW * swarm[i + n].getVelocity(j)
					+ lC1 * uniform<float>(rnd[n]) * ((*pbest[n])[j] - (*pos[n])[j]);

			rnd = rnds;

#pragma unroll
			for(size_t n = 0; n < 4; ++n)
			{
				v[n] = v[n] + lC2 * uniform<float>(rnd[n]) * (gbest[j] - (*pos[n])[j]);

				v[n] = clamp(v[n], vMax);

				swarm[i + n].setVelocity(j, v[n]);

				float S = 1.f / (1.f + expf(-v[n]));

				swarm[i + n].setPosition(j, rnds() < S ? 1 : 0);
			}
		}
	}
}

void MKPSolver::fitness(const Particle::arg_t * x[4], float * v)
{
#pragma unroll
	for(size_t p = 0; p < 4; ++p)
	{
		float f = 0;

		for(unsigned i = 0; i < data.objects(); ++i)
			f += data.profits()[i] * (*x[p])[i];
	
		float spl = 0;

		for(unsigned i = 0; i < data.knapsacks(); ++i)
		{
			float m = 0;

			for(unsigned j = 0; j < data.objects(); ++j)
				m += data.weights()[i * data.objects() + j] * (*x[p])[j];

			spl += poslin(penalty * (m - data.capacities()[i]));
		}

		f -= spl;

		v[p] = f;
	}
}

#ifdef SHOW_SYNC
	size_t MKPSolver::showStats(size_t n, bool sync)
#else
	size_t MKPSolver::showStats(size_t n)
#endif
{
	size_t worst = 0;

	mstats[0] = swarm[best].getBest();  // best
	mstats[1] = swarm[0].getBest();     // worst
	mstats[2] = swarm[0].getBest();     // average

	for(size_t i = 1; i < swarm.size(); ++i)
	{
		if(swarm[i].getBest() < mstats[1])
		{
			worst = i;
			mstats[1] = swarm[i].getBest();
		}

		mstats[2] += swarm[i].getBest();
	}

	mstats[2] /= swarm.size();

	MPI_Request request;

	if(id == MPIROOT)
	{
		cerr << endl << '[' << n << ']' << endl;
		cerr << "PROC\tBEST\t\tWORST\t\tAVERAGE" << endl;
		cerr << "===============================================" << endl;
		cerr << "0\t" << mstats[0] << "\t\t" << mstats[1] << "\t\t" << mstats[2] << endl;

		for(int i = 1; i < procs; ++i)
		{
			size_t offset = 0;

#ifdef SHOW_SYNC
			if(sync)
				MPI_Recv(mstats + offset, 3, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			else
#endif
			{
				offset = 3 * i;
				MPI_Irecv(mstats + offset, 3, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &request);
			}

			cerr << i << "\t";
			cerr << mstats[0 + offset] << "\t\t" << mstats[1 + offset] << "\t\t" << mstats[2 + offset] << endl;
		}
	}
	else
	{
#ifdef SHOW_SYNC
		if(sync)
			MPI_Send(mstats, 3, MPI_FLOAT, MPIROOT, 0, MPI_COMM_WORLD);
		else
#endif
			// mstats might contain incorrect data if overwritten before send,
			// but it doesn't matter.
			MPI_Isend(mstats, 3, MPI_FLOAT, MPIROOT, 0, MPI_COMM_WORLD, &request);
	}

	return worst;
}

bool MKPSolver::exchangeSync(size_t worst)
{
	struct
	{
		float val;
		int rank;
	} mpin, mpout;

	mpin.rank = id;

	mpin.val = (exLastFitness == swarm[best].getBest())
		? -std::numeric_limits<float>::infinity() : swarm[best].getBest();

	MPI_Allreduce(&mpin, &mpout, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);

	if(mpout.val == -std::numeric_limits<float>::infinity()) 	// probably converged
		return false;

	const Particle::arg_t & sol = swarm[best].solution();

	if(id == mpout.rank)
		memcpy(hsols, sol.data(), sol.bytes());

	MPI_Bcast(hsols, sol.bytes(), MPI_CHAR, mpout.rank, MPI_COMM_WORLD);

	if(id != mpout.rank)
	{
		swarm[worst].setBest(mpout.val, hsols);
		best = worst;
	}

	exLastFitness = mpout.val;

	return true;
}


MKPSolver::MKPSolver(const MKPData &dat, size_t particles)
:
	data(dat), swarm(adjustParticleCount(particles), data.objects()),
	dims(data.objects()), best(0)
{
	setInertiaWeight();
	setLearningFactors();
	setVelocityMax();
	setPenalty();

	init();
}

MKPSolver::~MKPSolver()
{
	delete [] hsols;
}

void MKPSolver::setInertiaWeight(float w)
{
	iW = w;
}

void MKPSolver::setLearningFactors(float c1, float c2)
{
	lC1 = c1;
	lC2 = c2;
}

void MKPSolver::setVelocityMax(float m)
{
	vMax = m;
}

void MKPSolver::setPenalty(float p)
{
	penalty = p;
}

const MKPSolver::sol_t & MKPSolver::solve(size_t epochs,
		size_t stop, size_t exchange, size_t statistics)
{
	assert(swarm.size() % 4 == 0);

	epochs = adjustEpochs(epochs);

	if(!statistics) statistics = epochs + 1;
	if(!exchange)   exchange   = epochs + 1;
	if(!stop)       stop       = epochs + 1;

	exLastFitness = std::numeric_limits<float>::infinity();

	size_t worst = 0;
	size_t m = stop;

#pragma omp parallel
	for(size_t n = 0, e = exchange - 1, s = statistics - 1; n < epochs; ++n, --e, --s)
	{
		stats.iterations = n + 1;

		findBest();

		if(!s)
		{
			s = statistics;

#pragma omp single
			worst = showStats(n);
		}

		if(!e)
		{
			e = exchange;

#pragma omp single
			{
				if(exchangeSync(worst)) m = stop;
				else                  --m;
			}
		}

		if(!m) break;

		update();
	}

	findBest();

	// Get best solution

	struct
	{
		float val;
		int rank;
	} mpin, mpout;

	mpin.rank = id;
	mpin.val = swarm[best].getBest();

	MPI_Allreduce(&mpin, &mpout, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);

	if(mpout.rank != MPIROOT)
	{
		const Particle::arg_t & sol = swarm[best].solution();

		if(id == mpout.rank)
		{
			memcpy(hsols, sol.data(), sol.bytes());
			MPI_Send(hsols, sol.bytes(), MPI_CHAR, MPIROOT, 1, MPI_COMM_WORLD);
		}
		else if(id == MPIROOT)
		{
			MPI_Recv(hsols, sol.bytes(), MPI_CHAR, mpout.rank, 1,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			swarm[worst].setBest(mpout.val, hsols);
			best = worst;
		}
	}

	return swarm[best].solution();
}

const MKPSolver::Stats & MKPSolver::statistics() const
{
	float worst = swarm[0].getBest();
	float avg = swarm[0].getBest();

	for(size_t i = 1; i < swarm.size(); ++i)
	{
		if(swarm[i].getBest() < worst)
			worst = swarm[i].getBest();

		avg += swarm[i].getBest();
	}

	stats.fitness.best = swarm[best].getBest();
	stats.fitness.worst = worst;
	stats.fitness.avg = avg / swarm.size();

	const Particle::arg_t &p = swarm[best].solution();

	for(size_t i = 0; i < data.knapsacks(); ++i)
	{
		size_t w = 0;

		for(size_t j = 0; j < dims; ++j)
			w += p[j] * data.weights()[i * data.objects() + j];

		stats.weight.push_back(w);
	}

	return stats;
}

