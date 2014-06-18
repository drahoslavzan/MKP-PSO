//===================================================================
// File:        knapsack.cu
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
#include "cuda.h"
#include "reduction.h"

#include <cassert>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <limits>

#include <mpi.h>


using namespace std;


// DEVICE

__device__ inline float clamp(float v, float m)
{
	if(v < -m) return -m;
	if(v > m) return m;
	return v;
}

__device__ inline float poslin(float v)
{
	return (v > 0.f) ? v : 0.f;
}


// KERNEL

KERNEL init(Swarm::Data swarm, float vMax, Seed4 seed)
{
	unsigned p = blockIdx.x * blockDim.x + threadIdx.x;

	if(p >= swarm.count()) return;

	Random4x32<float> rnd(p, 0xcaffe123, seed);

	for(unsigned d = 0; d < swarm.dims(); ++d)
	{
		swarm.setPosition(p, d, rnd() < 0.5f ? 1 : 0);
		swarm.velocity(p)[d] = rnd() * vMax;
	}
}

KERNEL update(Swarm::Data swarm, MKPSolver::ParamUpdate paramUpdate, Seed4 seed)
{
	unsigned p = blockIdx.x;
	unsigned d = threadIdx.x;

	assert(p < swarm.count());
	assert(d < swarm.dims());

	Random4x32<float> rnd(p * swarm.dims() + d, 0xabcf0198, seed);

	float v = paramUpdate.iW * swarm.velocity(p)[d]
		+ paramUpdate.lC1 * rnd() * (swarm.solution(p, d) - swarm.position(p, d))
		+ paramUpdate.lC2 * rnd() * (swarm.solution(swarm.best(), d)
				- swarm.position(p, d));

	v = clamp(v, paramUpdate.vMax);

	swarm.velocity(p)[d] = v;

	float S = 1.f / (1.f + expf(-v));
 	
	swarm.setPosition(p, d, (rnd() < S) ? 1 : 0);
}

template <unsigned blockSize>
KERNEL fitnessShared(Swarm::Position x, MKPSolver::ParamData data)
{
	assert(data.n == blockSize);

	unsigned p = blockIdx.x;
	unsigned d = threadIdx.x;

	assert(p < x.count());
	assert(d < blockSize);

	__shared__ float sdata[blockSize];

	sdata[d] = data.d_pro[d] * x.position(p, d);

	rdSum<blockSize>(sdata);

	__shared__ float f[2];

	if(d == 0)
	{
		f[0] = sdata[0];
		f[1] = 0.f;
	}

	for(unsigned i = 0; i < data.m; ++i)
	{
		sdata[d] = data.d_wei[i * blockSize + d] * x.position(p, d);

		rdSum<blockSize>(sdata);

		if(d == 0)
			f[1] += poslin(data.penalty * (sdata[0] - data.d_cap[i]));
	}

	if(d == 0)
		x.setBest(p, f[0] - f[1]);
}


// HOST

size_t MKPSolver::adjustParticleCount(size_t n)
{
	const size_t N = PSO_MIN_PARTICLES;

	return (n % N) ? n + (N - (n % N)) : n;
}

size_t MKPSolver::adjustEpochs(size_t n)
{
	return n;
}

#ifdef SHOW_SYNC
	void MKPSolver::showStats(size_t n, bool sync)
#else
	void MKPSolver::showStats(size_t n)
#endif
{
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
			MPI_Send(mstats, 3, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
		else
#endif
			MPI_Isend(mstats, 3, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &request);
	}
}

void MKPSolver::exchangeFinal(float best)
{
	struct
	{
		float val;
		int rank;
	} mpin, mpout;

	mpin.rank = id;
	mpin.val = best;

	MPI_Allreduce(&mpin, &mpout, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);

	if(mpout.rank != MPIROOT)
	{
		if(id == mpout.rank)
		{
			swarm.getBest(h_sol);
			MPI_Send(h_sol, swarm.data().dims(), MPI_CHAR /* bit_t */, MPIROOT, 1, MPI_COMM_WORLD);
		}
		else if(id == MPIROOT)
		{
			MPI_Recv(h_sol, swarm.data().dims(), MPI_CHAR, mpout.rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			swarm.setBest(mpout.val, h_sol);
		}
	}
}


MKPSolver::MKPSolver(const MKPData &dat, size_t particles)
:
	mkpdata(dat), swarm(adjustParticleCount(particles), dat.objects(), dat.align()),
	h_sol(NULL)
{
	paramData.d_pro = NULL;
	paramData.d_wei = NULL;
	paramData.d_cap = NULL;

	setInertiaWeight();
	setLearningFactors();
	setVelocityMax();
	setPenalty();

	init();
}

MKPSolver::~MKPSolver()
{
	cudaFree(paramData.d_pro);
	cudaFree(paramData.d_wei);
	cudaFree(paramData.d_cap);

	delete [] h_sol;
}

void MKPSolver::setInertiaWeight(float w)
{
	paramUpdate.iW = w;
}

void MKPSolver::setLearningFactors(float c1, float c2)
{
	paramUpdate.lC1 = c1;
	paramUpdate.lC2 = c2;
}

void MKPSolver::setVelocityMax(float m)
{
	paramUpdate.vMax = m;
}

void MKPSolver::setPenalty(float p)
{
	paramData.penalty = p;
}

void MKPSolver::init()
{
	size_t size = sizeof(float);

	paramData.n = mkpdata.align();
	paramData.m = mkpdata.knapsacks();

	checkCudaErrors(cudaMalloc(&paramData.d_pro, paramData.n * size));
	checkCudaErrors(cudaMalloc(&paramData.d_wei,
				paramData.n * paramData.m * size));
	checkCudaErrors(cudaMalloc(&paramData.d_cap, paramData.m * size));

	cudaMemcpy(paramData.d_pro, mkpdata.profits(), paramData.n * size,
			cudaMemcpyHostToDevice);
	cudaMemcpy(paramData.d_wei, mkpdata.weights(),
			paramData.n * paramData.m * size, cudaMemcpyHostToDevice);
	cudaMemcpy(paramData.d_cap, mkpdata.capacities(), paramData.m * size,
			cudaMemcpyHostToDevice);

	DEVICE(::init, swarm.data().count())(swarm.data(), paramUpdate.vMax, Seed4());

	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	h_sol = new char[swarm.data().dims()];

	size_t mstatssz;

	if(id == MPIROOT) mstatssz = 3 * procs;
	else mstatssz = 3;

	mstats = new float[mstatssz];
	memset(mstats, 0, mstatssz);
}

template <unsigned blockSize>
Swarm::sol_t MKPSolver::solveSync(size_t epochs, size_t stop,
		size_t exchange, size_t statistics)
{
	Swarm::Data &data = swarm.data();

	assert(paramData.n == blockSize);
	assert(paramData.n == mkpdata.align());
	assert(data.dims() == mkpdata.objects());

	exLastFitness = std::numeric_limits<float>::infinity();

	struct
	{
		float val;
		int rank;
	} mpin, mpout;

	mpin.rank = id;

	size_t worst = 0;

	for(size_t n = 0, dummy, m = stop, e = exchange - 1, s = statistics - 1;
			n < epochs; ++n, --e, --s)
	{
		stats.iterations = n + 1;

		fitnessShared<blockSize><<<data.count(), blockSize>>>(data.position(), paramData);

		if(!s)
		{
			s = statistics;

			swarm.findBest(mstats, &worst, &dummy);

			mpin.val = mstats[0];

			showStats(n);
		}
		else
			mpin.val = swarm.findBest(&dummy);

		if(!e)
		{
			e = exchange;

			float f = mpin.val;

			if(mpin.val == exLastFitness)
				mpin.val = -std::numeric_limits<float>::infinity();

			MPI_Allreduce(&mpin, &mpout, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);

			if(mpout.val == -std::numeric_limits<float>::infinity())
			{
				--m;
				mpout.val = f;
			}
			else
			{
				m = stop;
				exLastFitness = mpout.val;

				if(id == mpout.rank)
					swarm.getBest(h_sol);
	
				MPI_Bcast(h_sol, data.dims(), MPI_CHAR /* bit_t */, mpout.rank,
						MPI_COMM_WORLD);

				if(id != mpout.rank)
					swarm.setBest(mpout.val, h_sol, worst);
			}
		}
		else
			mpout.val = mpin.val;

		if(!m) break;

		update<<<data.count(), data.dims()>>>(data, paramUpdate, Seed4());
	}

	exchangeFinal(mpout.val);

	return swarm.solution();
}

Swarm::sol_t MKPSolver::solve(size_t epochs, size_t stop, size_t exchange,
		size_t statistics)
{
	assert(paramData.n == swarm.data().stride());

	epochs = adjustEpochs(epochs);

	if(!stop)       stop       = epochs + 1;
	if(!exchange)   exchange   = epochs + 1;
	if(!statistics) statistics = epochs + 1;

	switch(paramData.n)
	{
//#if __CUDA_ARCH__ >= 200
		case 1024:
			return solveSync<1024>(epochs, stop, exchange, statistics);
//#endif
		case 512:
			return solveSync< 512>(epochs, stop, exchange, statistics);
		case 256:
			return solveSync< 256>(epochs, stop, exchange, statistics);
		case 128:
			return solveSync< 128>(epochs, stop, exchange, statistics);
		case 64:
			return solveSync<  64>(epochs, stop, exchange, statistics);
		case 32:
			return solveSync<  32>(epochs, stop, exchange, statistics);
		case 16:
			return solveSync<  16>(epochs, stop, exchange, statistics);
		case 8:
			return solveSync<   8>(epochs, stop, exchange, statistics);
		case 4:
			return solveSync<   4>(epochs, stop, exchange, statistics);
		case 2:
			return solveSync<   2>(epochs, stop, exchange, statistics);
		case 1:
			return solveSync<   1>(epochs, stop, exchange, statistics);
		default:
		{
			stringstream iss;
			iss << paramData.n;
			throw std::string("solve(): ") + iss.str() + ": Invalid block size";
		}
	}
}

const MKPSolver::Stats & MKPSolver::statistics() const
{
	Swarm::Host host(swarm.data().count(), swarm.data().dims(), swarm.data().stride());

	host << swarm;

	Swarm::Position &pos = host.data().position();

	float worst = pos.best(0);
	float avg = worst;

	for(size_t i = 1; i < pos.count(); ++i)
	{
		if(pos.best(i) < worst)
			worst = pos.best(i);

		avg += pos.best(i);
	}

	stats.fitness.best = pos.best(host.data().best());
	stats.fitness.worst = worst;
	stats.fitness.avg = avg / pos.count();

	size_t hb = host.data().best();

	for(size_t i = 0; i < mkpdata.knapsacks(); ++i)
	{
		size_t w = 0;

		for(size_t j = 0; j < pos.dims(); ++j)
			w += pos.solution(hb, j) * mkpdata.weights()[i * mkpdata.align() + j];

		stats.weight.push_back(w);
	}

	return stats;
}

