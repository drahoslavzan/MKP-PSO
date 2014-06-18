//===================================================================
// File:        reduction.h
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    Reduction kernels.
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


#ifndef _REDUCTION_H_
#define _REDUCTION_H_


#ifdef __CUDACC__
#  include "cuda.h"
#endif


struct ValPos
{
	float val;
	unsigned pos;
};


#ifdef __CUDACC__
	template <unsigned blockSize>
	__device__ void rdSum(float * sdata)
	{
		unsigned tid = threadIdx.x;

		__syncthreads();

#  if __CUDA_ARCH__ >= 200
		if(blockSize >= 1024)
		{
			if(tid < 512)
				sdata[tid] += sdata[tid + 512];
			__syncthreads();
		}
#  endif
		if(blockSize >= 512)
		{
			if(tid < 256)
				sdata[tid] += sdata[tid + 256];
			__syncthreads();
		}
		if(blockSize >= 256)
		{
			if(tid < 128)
				sdata[tid] += sdata[tid + 128];
			__syncthreads();
		}
		if(blockSize >= 128)
		{
			if(tid < 64)
				sdata[tid] += sdata[tid +  64];
			__syncthreads();
		}

		if(tid < 32)
		{
			volatile float *smem = sdata;

			if(blockSize >= 64) { smem[tid] += smem[tid + 32]; }
			if(blockSize >= 32) { smem[tid] += smem[tid + 16]; }
			if(blockSize >= 16) { smem[tid] += smem[tid +  8]; }
			if(blockSize >=  8) { smem[tid] += smem[tid +  4]; }
			if(blockSize >=  4) { smem[tid] += smem[tid +  2]; }
			if(blockSize >=  2) { smem[tid] += smem[tid +  1]; }
		}
	}
#endif

#ifdef __CUDACC__
	template <unsigned blockSize>
	KERNEL rdMax(float * vec, unsigned sz, ValPos * max, unsigned leader)
	{
		__shared__ float sdata[blockSize];

		unsigned tid = threadIdx.x;

		if(tid < sz)
			sdata[tid] = vec[tid];
		else
			sdata[tid] = -FLT_MAX;

		__syncthreads();

		float oldVal;

		if(tid == 0) oldVal = sdata[leader];

#  if __CUDA_ARCH__ >= 200
		if(blockSize >= 1024)
		{
			if(tid < 512)
				sdata[tid] = fmaxf(sdata[tid], sdata[tid + 512]);
			__syncthreads();
		}
#  endif
		if(blockSize >= 512)
		{
			if(tid < 256)
				sdata[tid] = fmaxf(sdata[tid], sdata[tid + 256]);
			__syncthreads();
		}
		if(blockSize >= 256)
		{
			if(tid < 128)
				sdata[tid] = fmaxf(sdata[tid], sdata[tid + 128]);
			__syncthreads();
		}
		if(blockSize >= 128)
		{
			if(tid < 64)
				sdata[tid] = fmaxf(sdata[tid], sdata[tid +  64]);
			__syncthreads();
		}

		if(tid < 32)
		{
			volatile float *smem = sdata;

			if(blockSize >= 64) { smem[tid] = fmaxf(smem[tid], smem[tid + 32]); }
			if(blockSize >= 32) { smem[tid] = fmaxf(smem[tid], smem[tid + 16]); }
			if(blockSize >= 16) { smem[tid] = fmaxf(smem[tid], smem[tid +  8]); }
			if(blockSize >=  8) { smem[tid] = fmaxf(smem[tid], smem[tid +  4]); }
			if(blockSize >=  4) { smem[tid] = fmaxf(smem[tid], smem[tid +  2]); }
			if(blockSize >=  2) { smem[tid] = fmaxf(smem[tid], smem[tid +  1]); }
		}

		__syncthreads();

		// WARN: might be dangerous
		if(tid < sz && vec[tid] == sdata[0]) max[0].pos = tid;

		if(tid == 0)
		{
			max[0].val = sdata[0];

			// new best
			max[1].pos = oldVal != sdata[0];
		}
	}

	template <unsigned blockSize>
	KERNEL rdMaxMinSum(float * vec, unsigned sz, ValPos * stats, unsigned leader)
	{
		__shared__ float sdata[3 * blockSize];

		float * smax = sdata + 0 * blockSize;
		float * smin = sdata + 1 * blockSize;
		float * ssum = sdata + 2 * blockSize;

		unsigned tid = threadIdx.x;

		if(tid < sz)
		{
			smax[tid] = vec[tid];
			smin[tid] = vec[tid];
			ssum[tid] = vec[tid];
		}
		else
		{
			smax[tid] = -FLT_MAX;
			smin[tid] =  FLT_MAX;
			ssum[tid] =  0;
		}

		__syncthreads();

		float oldVal;

		if(tid == 0) oldVal = smax[leader];

#  if __CUDA_ARCH__ >= 200
		if(blockSize >= 1024)
		{
			if(tid < 512)
			{
				smax[tid]  = fmaxf(smax[tid], smax[tid + 512]);
				smin[tid]  = fminf(smin[tid], smin[tid + 512]);
				ssum[tid] += ssum[tid + 512];
			}
			__syncthreads();
		}
#  endif
		if(blockSize >= 512)
		{
			if(tid < 256)
			{
				smax[tid]  = fmaxf(smax[tid], smax[tid + 256]);
				smin[tid]  = fminf(smin[tid], smin[tid + 256]);
				ssum[tid] += ssum[tid + 256];
			}
			__syncthreads();
		}
		if(blockSize >= 256)
		{
			if(tid < 128)
			{
				smax[tid]  = fmaxf(smax[tid], smax[tid + 128]);
				smin[tid]  = fminf(smin[tid], smin[tid + 128]);
				ssum[tid] += ssum[tid + 128];
			}
			__syncthreads();
		}
		if(blockSize >= 128)
		{
			if(tid < 64)
			{
				smax[tid]  = fmaxf(smax[tid], smax[tid +  64]);
				smin[tid]  = fminf(smin[tid], smin[tid +  64]);
				ssum[tid] += ssum[tid + 64];
			}
			__syncthreads();
		}

		if(tid < 32)
		{
			volatile float *vmax = smax;
			volatile float *vmin = smin;
			volatile float *vsum = ssum;

			if(blockSize >= 64)
			{
				vmax[tid]  = fmaxf(vmax[tid], vmax[tid + 32]);
				vmin[tid]  = fminf(vmin[tid], vmin[tid + 32]);
				vsum[tid] += vsum[tid + 32];
			}
			if(blockSize >= 32)
			{
				vmax[tid]  = fmaxf(vmax[tid], vmax[tid + 16]);
				vmin[tid]  = fminf(vmin[tid], vmin[tid + 16]);
				vsum[tid] += vsum[tid + 16];
			}
			if(blockSize >= 16)
			{
				vmax[tid]  = fmaxf(vmax[tid], vmax[tid +  8]);
				vmin[tid]  = fminf(vmin[tid], vmin[tid +  8]);
				vsum[tid] += vsum[tid + 8];
			}
			if(blockSize >=  8)
			{
				vmax[tid]  = fmaxf(vmax[tid], vmax[tid +  4]);
				vmin[tid]  = fminf(vmin[tid], vmin[tid +  4]);
				vsum[tid] += vsum[tid + 4];
			}
			if(blockSize >=  4)
			{
				vmax[tid]  = fmaxf(vmax[tid], vmax[tid +  2]);
				vmin[tid]  = fminf(vmin[tid], vmin[tid +  2]);
				vsum[tid] += vsum[tid + 2];
			}
			if(blockSize >=  2)
			{
				vmax[tid]  = fmaxf(vmax[tid], vmax[tid +  1]);
				vmin[tid]  = fminf(vmin[tid], vmin[tid +  1]);
				vsum[tid] += vsum[tid + 1];
			}
		}

		__syncthreads();

		// WARN: might be dangerous
		if(tid < sz)
		{
			if(vec[tid] == smax[0]) stats[0].pos = tid;
			if(vec[tid] == smin[0]) stats[1].pos = tid;
		}

		if(tid == 0)
		{
			stats[0].val = smax[0];
			stats[1].val = smin[0];
			stats[2].val = ssum[0];

			// new best
			stats[2].pos = oldVal != smax[0];
		}
	}

	static KERNEL fillt(float * vec, float val)
	{
		unsigned tid = threadIdx.x;

		vec[tid] = val;
	}
#endif


#endif /* _REDUCTION_H_ */
