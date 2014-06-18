//===================================================================
// File:        swarm.h
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    The PSO swarm.
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


#ifndef _SWARM_H_
#define _SWARM_H_


#include "cuda.h"
#include "reduction.h"
#include "generic.h"

#include <limits>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <cmath>


namespace Swarm
{
	class Host;
	class Device;
	class Position;

	const Host & operator <<(Device &d, const Host &h);
	const Device & operator <<(Host &h, const Device &d);

	typedef char bit_t;
	typedef std::vector<bit_t> sol_t;

	class Position
	{
		friend class Host;
		friend class Device;
		friend class Data;
		friend const Host & operator <<(Device &d, const Host &h);
		friend const Device & operator <<(Host &h, const Device &d);
	private:
		const size_t cnt, dms, str;
		float * bst;
		bit_t * pos, * sol;
	private:
#ifndef __CUDACC__
		Position(size_t count, size_t dims, size_t align);
#else
		__host__ Position(size_t count, size_t dims, size_t align)
		:
			cnt(count), dms(dims), str(align), bst(NULL), pos(NULL), sol(NULL)
		{
			assert(align > dims);
		}
	public:
		__host__ __device__ inline bit_t position(size_t p, size_t d)
		{
			assert(p < cnt);

			return (pos + p * str)[d];
		}

		__host__ __device__ inline void setPosition(size_t p, size_t d, bit_t v)
		{
			assert(p < cnt);

			(pos + p * str)[d] = v;
		}

		__host__ __device__ inline bit_t solution(size_t p, size_t d) const
		{
			assert(p < cnt);

			return (sol + p * str)[d];
		}

		__host__ __device__ inline float best(size_t i) const
		{
			assert(i < cnt);

			return bst[i];
		}

		__host__ __device__ inline void setBest(size_t i, float v)
		{
			assert(i < cnt);

			if(v > bst[i])
			{
				bst[i] = v;

				memcpy(sol + i * str, pos + i * str, dms * sizeof(bit_t));
			}
		}
	public:
		__host__ __device__ inline size_t count() const
		{
			return cnt;
		}
		__host__ __device__ inline size_t dims() const
		{
			return dms;
		}
		__host__ __device__ inline size_t stride() const
		{
			return str;
		}
		/*
		__host__ __device__ inline size_t size() const
		{
			return cnt * dms;
		}
		*/
#endif
	};

	class Data
	{
		friend class Host;
		friend class Device;
		friend const Host & operator <<(Device &d, const Host &h);
		friend const Device & operator <<(Host &h, const Device &d);
	private:
		size_t gbest;
		float * vel;
		Position pos;
	private:
#ifndef __CUDACC__
		Data(size_t count, size_t dims, size_t align);
#else
		__host__ Data(size_t count, size_t dims, size_t align)
		:
			gbest(0), pos(count, dims, align)
		{
			assert(align > dims);
		}
	public:
		__host__ __device__ inline float * velocity(size_t i)
		{
			assert(i < count());

			return vel + i * stride();
		}

		__host__ __device__ inline bit_t solution(size_t p, size_t d) const
		{
			assert(p < count());

			return pos.solution(p, d);
		}

		__host__ __device__ inline bit_t position(size_t p, size_t d)
		{
			assert(p < count());

			return pos.position(p, d);
		}

		__host__ __device__ inline void setPosition(size_t p, size_t d, bit_t v)
		{
			assert(p < count());

			pos.setPosition(p, d, v);
		}
	public:
		__host__ __device__ inline size_t best() const
		{
			return gbest;
		}
		__host__ __device__ inline size_t count() const
		{
			return pos.count();
		}
		__host__ __device__ inline size_t dims() const
		{
			return pos.dims();
		}
		__host__ __device__ inline size_t stride() const
		{
			return pos.stride();
		}
		/*
		__host__ __device__ inline size_t size() const
		{
			return pos.size();
		}
		*/
	public:
		__host__ inline Position & position()
		{
			return pos;
		}
#endif
	};

	class Device
	{
		friend const Host & operator <<(Device &d, const Host &h);
		friend const Device & operator <<(Host &h, const Device &d);
	private:
		Data d_data;
		ValPos * d_stats;
		size_t np2cnt;
	public:
#ifndef __CUDACC__
		Device(size_t count, size_t dims, size_t align);
#else
		__host__ Device(size_t count, size_t dims, size_t align)
		:
			d_data(count, dims, align), np2cnt(nearestPow2(count))
		{
			assert(align > dims);

			size_t size = count * align * sizeof(float);
			size_t sizeB = count * align * sizeof(bit_t);

			checkCudaErrors(cudaMalloc(&d_stats, 3 * sizeof(ValPos)));

			checkCudaErrors(cudaMalloc(&d_data.vel, size));
			checkCudaErrors(cudaMalloc(&d_data.pos.sol, sizeB));
			checkCudaErrors(cudaMalloc(&d_data.pos.pos, sizeB));

			checkCudaErrors(cudaMalloc(&d_data.pos.bst, count * sizeof(float)));
			fillt<<<1, count>>>(d_data.pos.bst, -std::numeric_limits<float>::max());
		}

		__host__ ~Device()
		{
			cudaFree(d_data.vel);
			cudaFree(d_data.pos.sol);
			cudaFree(d_data.pos.pos);
			cudaFree(d_data.pos.bst);
			cudaFree(d_stats);
		}

		__host__ inline Data & data()
		{
			return d_data;
		}

		__host__ inline const Data & data() const
		{
			return d_data;
		}

		__host__ inline void getBest(bit_t * best)
		{
			cudaMemcpy(best, d_data.pos.sol + d_data.best() * d_data.stride(),
					d_data.dims() * sizeof(bit_t), cudaMemcpyDeviceToHost);
		}

		__host__ inline void setBest(float fit, bit_t * best, size_t i = 0)
		{
			assert(i < d_data.count());

			d_data.gbest = i;

			cudaMemcpy(d_data.pos.bst + i, &fit, sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(d_data.pos.pos + i * d_data.stride(), best,
					d_data.dims() * sizeof(bit_t), cudaMemcpyHostToDevice);
			cudaMemcpy(d_data.pos.sol + i * d_data.stride(),
					d_data.pos.pos + i * d_data.stride(),
					d_data.dims() * sizeof(bit_t), cudaMemcpyDeviceToDevice);
		}

		__host__ inline void findBest(float * sts, size_t * worst, size_t * nlead)
		{
			switch(np2cnt)
			{
//#  if __CUDA_ARCH__ >= 200
				case 1024:
					rdMaxMinSum<1024><<<1, 1024>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
//#  endif
				case 512:
					rdMaxMinSum< 512><<<1,  512>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 256:
					rdMaxMinSum< 256><<<1,  256>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 128:
					rdMaxMinSum< 128><<<1,  128>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 64:
					rdMaxMinSum<  64><<<1,   64>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 32:
					rdMaxMinSum<  32><<<1,   32>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 16:
					rdMaxMinSum<  16><<<1,   16>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 8:
					rdMaxMinSum<   8><<<1,    8>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 4:
					rdMaxMinSum<   4><<<1,    4>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 2:
					rdMaxMinSum<   2><<<1,    2>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
			}

			ValPos stats[3];

			cudaMemcpy(stats, d_stats, 3 * sizeof(ValPos), cudaMemcpyDeviceToHost);

			sts[0] = stats[0].val;
			sts[1] = stats[1].val;
			sts[2] = stats[2].val / d_data.count();

			d_data.gbest = stats[0].pos;
			*worst = stats[1].pos;

			if(stats[2].pos) *nlead = 0;
		}

		__host__ inline float findBest(size_t * nlead)
		{
			switch(np2cnt)
			{
//#  if __CUDA_ARCH__ >= 200
				case 1024:
					rdMax<1024><<<1, 1024>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
//#  endif
				case 512:
					rdMax< 512><<<1,  512>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 256:
					rdMax< 256><<<1,  256>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 128:
					rdMax< 128><<<1,  128>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 64:
					rdMax<  64><<<1,   64>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 32:
					rdMax<  32><<<1,   32>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 16:
					rdMax<  16><<<1,   16>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 8:
					rdMax<   8><<<1,    8>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 4:
					rdMax<   4><<<1,    4>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
				case 2:
					rdMax<   2><<<1,    2>>>(d_data.pos.bst, d_data.count(), d_stats, d_data.gbest);
					break;
			}

			ValPos max[2];

			cudaMemcpy(max, d_stats, 2 * sizeof(ValPos), cudaMemcpyDeviceToHost);

			d_data.gbest = max[0].pos;

			if(max[1].pos) *nlead = 0;

			return max[0].val;
		}

		__host__ inline sol_t solution() const
		{
			sol_t sol(d_data.dims());

			cudaMemcpy(sol.data(), d_data.pos.sol + d_data.best() * d_data.stride(),
					d_data.dims() * sizeof(bit_t), cudaMemcpyDeviceToHost);

			return sol;
		}
#endif
	};

	class Host
	{
		friend const Host & operator <<(Device &d, const Host &h);
		friend const Device & operator <<(Host &h, const Device &d);
	private:
		Data h_data;
	public:
#ifndef __CUDACC__
		Host(size_t count, size_t dims, size_t align);
#else
		__host__ Host(size_t count, size_t dims, size_t align)
		:
			h_data(count, dims, align)
		{
			assert(align > dims);

			size_t size = count * align;

			h_data.vel = new float[size];
			h_data.pos.sol = new bit_t[size];
			h_data.pos.pos = new bit_t[size];

			h_data.pos.bst = new float[count];

			for(size_t i = 0; i < count; ++i)
				h_data.pos.bst[i] = -std::numeric_limits<float>::max();
		}
		__host__ ~Host()
		{
			delete [] h_data.vel;
			delete [] h_data.pos.sol;
			delete [] h_data.pos.pos;
			delete [] h_data.pos.bst;
		}

		__host__ inline Data & data()
		{
			return h_data;
		}

		__host__ inline const Data & data() const
		{
			return h_data;
		}

		__host__ inline sol_t solution() const
		{
			sol_t sol(h_data.dims());

			memcpy(sol.data(), h_data.pos.sol + h_data.best() * h_data.stride(),
					h_data.dims() * sizeof(bit_t));

			return sol;
		}
#endif
	};

#ifdef __CUDACC__
	__host__ inline const Host & operator <<(Device &d, const Host &h)
	{
		size_t size = h.h_data.count() * h.h_data.stride() * sizeof(float);
		size_t sizeB = h.h_data.count() * h.h_data.stride() * sizeof(bit_t);

		d.d_data.gbest = h.h_data.gbest;

		cudaMemcpy(d.d_data.vel, h.h_data.vel, size, cudaMemcpyHostToDevice);
		cudaMemcpy(d.d_data.pos.sol, h.h_data.pos.sol, sizeB, cudaMemcpyHostToDevice);
		cudaMemcpy(d.d_data.pos.pos, h.h_data.pos.pos, sizeB, cudaMemcpyHostToDevice);
		cudaMemcpy(d.d_data.pos.bst, h.h_data.pos.bst,
				h.h_data.count() * sizeof(float), cudaMemcpyHostToDevice);

		return h;
	}

	__host__ inline const Device & operator <<(Host &h, const Device &d)
	{
		size_t size = d.d_data.count() * d.d_data.stride() * sizeof(float);
		size_t sizeB = h.h_data.count() * h.h_data.stride() * sizeof(bit_t);

		h.h_data.gbest = d.d_data.gbest;

		cudaMemcpy(h.h_data.vel, d.d_data.vel, size, cudaMemcpyDeviceToHost);
		cudaMemcpy(h.h_data.pos.sol, d.d_data.pos.sol, sizeB, cudaMemcpyDeviceToHost);
		cudaMemcpy(h.h_data.pos.pos, d.d_data.pos.pos, sizeB, cudaMemcpyDeviceToHost);
		cudaMemcpy(h.h_data.pos.bst, d.d_data.pos.bst,
				d.d_data.count() * sizeof(float), cudaMemcpyDeviceToHost);

		return d;
	}
#endif
}


#endif /* _SWARM_H_ */
