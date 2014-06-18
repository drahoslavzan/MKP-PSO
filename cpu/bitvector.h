//===================================================================
// File:        bitvector.h
// Author:      Drahoslav Zan
// Email:       izan@fit.vutbr.cz
// Affiliation: Brno University of Technology,
//              Faculty of Information Technology
// Date:        Tue Nov 17 22:12:37 CET 2013
// Comments:    Container for binary values.
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


#ifndef _BITVECTOR_H_
#define _BITVECTOR_H_


#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <string>


template <typename T = unsigned char>
class bitvector
{
public:

	typedef T elem_t;
	static const size_t ELEM_BITS = 8 * sizeof(elem_t);

private:

	elem_t * array;
	size_t bitlen, elemlen, alloclen;

private:

	inline size_t bitLen2elemLen(size_t n)
	{
		return (n + (ELEM_BITS - 1)) / ELEM_BITS;
	}

	inline void alloc(size_t n)
	{
		if(n <= alloclen)
			return;

		free();
		array = (elem_t *) malloc(n * sizeof(elem_t));
		alloclen = n;
	}

	inline void realloc(size_t n)
	{
		if(n <= alloclen)
			return;

		array = (elem_t *) ::realloc(array, n * sizeof(elem_t));
		alloclen = n;
	}

	inline void free()
	{
		::free(array);
	}

public:

	bitvector(size_t n = 0)
		: array(NULL), bitlen(n), elemlen(bitLen2elemLen(n)), alloclen(0)
	{
		alloc(elemlen);
	}

	bitvector(size_t n, int v)
		: array(NULL), bitlen(n), elemlen(bitLen2elemLen(n)), alloclen(0)
	{
		assert(elemlen >= 1);
		assert(v == 0 || v == 1);

		alloc(elemlen);
		memset(array, v * 0xFF, sizeof(elem_t) * elemlen);
	}

	bitvector(const bitvector &bv)
		: array(NULL), alloclen(0)
	{
		operator =(bv);
	}

	~bitvector()
	{
		free();
	}

public:

	inline void assign(const elem_t *d)
	{
		operator =(d);
	}

	void assign(const elem_t *d, size_t n)
	{
		bitlen = n;
		elemlen = bitLen2elemLen(n);

		alloc(elemlen);

		assign(d);
	}

	void resize(size_t s, int p = 0)
	{
		assert(s > 0);
		assert(p == 0 || p == 1);

		size_t b = bitLen2elemLen(s);

		realloc(b);

		if(s > bitlen)
		{
			if(b > elemlen)
				memset(array + elemlen, p * 0xFF, (b - elemlen) * sizeof(elem_t));

			size_t m = bitlen / ELEM_BITS;
			size_t n = bitlen % ELEM_BITS;

			array[m] = (array[m] & ~(~0 << n)) | (p * (~0 << n));
		}

		bitlen = s;
		elemlen = b;
	}

	int get(size_t i) const
	{
		assert(i < bitlen);

		size_t b = i / ELEM_BITS;
		size_t s = i % ELEM_BITS;

		return (array[b] & (1 << s)) >> s;
	}

	void set(size_t i, int v = 1)
	{
		assert(i < bitlen);
		assert(v == 0 || v == 1);

		size_t b = i / ELEM_BITS;
		size_t s = i % ELEM_BITS;

		array[b] = (array[b] & ~(1 << s)) | (v << s);
	}

	void flip(size_t i)
	{
		assert(i < bitlen);

		size_t b = i / ELEM_BITS;
		size_t s = i % ELEM_BITS;

		array[b] = array[b] ^ (1 << s);
	}

	void push_back(int v)
	{
		elemlen = bitLen2elemLen(++bitlen);

		realloc(elemlen);

		set(bitlen - 1, v);
	}

	int pop_back()
	{
		assert(bitlen > 0);

		int v = get(bitlen - 1);

		elemlen = bitLen2elemLen(--bitlen);

		return v;
	}

	inline size_t size() const { return bitlen; }
	inline size_t length() const { return size(); }
	inline size_t bits() const { return size(); }

	inline size_t elems() const { return elemlen; }
	inline size_t bytes() const { return sizeof(elem_t) * elems(); }

	inline elem_t * data() { return array; }
	inline const elem_t * data() const { return array; }

	inline int operator[](size_t i) const { return get(i); }
	inline int operator()(size_t i) const { return get(i); }
	inline void operator()(size_t i, int v) { set(i, v); }

	inline const elem_t * operator =(const elem_t * d)
	{
		memcpy(array, d, sizeof(elem_t) * elemlen);
		
		return d;
	}

	const bitvector & operator =(const bitvector & bv)
	{
		bitlen = bv.bitlen;
		elemlen = bv.elemlen;

		alloc(bv.elemlen);

		memcpy(array, bv.array, sizeof(elem_t) * bv.elemlen);
		
		return bv;
	}

	inline operator elem_t * () { return data(); }
	inline operator const elem_t * () const { return data(); }

	operator std::string() const
	{
		std::string str;

		if(!size())
			return str;

		str.reserve(size());

		for(size_t i = 0; i < size(); ++i)
			str += '0' + get(i);
	
		return str;
	}

public:

	friend std::ostream & operator <<(std::ostream & os, const bitvector & bv)
	{
		if(!bv.size())
			return os;
	
		os << bv[0];
	
		for(size_t i = 1; i < bv.size(); ++i)
		{
			if(!(i % 64)) os << std::endl;
			else if(!(i % bitvector::ELEM_BITS)) os << " ";
	
			os << bv[i];
		}
	
		return os;
	}
};

typedef bitvector<> BitVector;


#endif /* _BITVECTOR_H_ */
