#pragma once

#include "cuda_utils.cuh"
#include "../../lisflood.h"
#include "index_1D.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct AssembledSolution
{
	NUMERIC_TYPE*     h;
//	NUMERIC_TYPE*     qx;
//	NUMERIC_TYPE*     qy;
//	NUMERIC_TYPE*     q3;
//	NUMERIC_TYPE*     q4;
	NUMERIC_TYPE*     z0;
	NUMERIC_TYPE*     z1x;
	NUMERIC_TYPE*     z1y;
//	NUMERIC_TYPE*     zxy;
	NUMERIC_TYPE*     n0;
	NUMERIC_TYPE* v;

	index_1D* act_idcs;
	int*      levels;
	int*      nghbr_counts;
	int*      cumu_nghbr_counts;
	int       length;
	bool      is_copy = false;

	AssembledSolution(const int& num_finest_elems, bool non_uniform_n)
	{
		size_t bytes_real = num_finest_elems * sizeof(NUMERIC_TYPE);
		size_t bytes_int  = num_finest_elems * sizeof(index_1D);

		h                 = (NUMERIC_TYPE*)malloc_device(bytes_real);
		v                 = (NUMERIC_TYPE*)malloc_device(bytes_real);
//		qx                = (NUMERIC_TYPE*)malloc_device(bytes_real);
//		qy                = (NUMERIC_TYPE*)malloc_device(bytes_real);
//		q3                = (NUMERIC_TYPE*)malloc_device(bytes_real);
//		q4                = (NUMERIC_TYPE*)malloc_device(bytes_real);
		z0                = (NUMERIC_TYPE*)malloc_device(bytes_real);
		z1x               = (NUMERIC_TYPE*)malloc_device(bytes_real);
		z1y               = (NUMERIC_TYPE*)malloc_device(bytes_real);
//		zxy               = (NUMERIC_TYPE*)malloc_device(bytes_real);


		n0                = (non_uniform_n) ? (NUMERIC_TYPE*)malloc_device(bytes_real) : nullptr;

		act_idcs          = (index_1D*)malloc_device(bytes_int);
		levels            = (int*)malloc_device(bytes_int);
		nghbr_counts      = (int*)malloc_device(bytes_int);
		cumu_nghbr_counts = (int*)malloc_device(bytes_int);
		length            = num_finest_elems;
//		is_copy           = false;
	}

	AssembledSolution(const AssembledSolution& original) { *this = original; is_copy = true; }

	~AssembledSolution()
	{
		if (!is_copy)
		{
			free_device(h);
//			free_device(qx);
//			free_device(qy);
//			free_device(q3);
//			free_device(q4);
			free_device(z0);
			free_device(z1x);
			free_device(z1y);
			free_device(v);

//			free_device(zxy);
			free_device(n0);
			free_device(act_idcs);
			free_device(levels);
			free_device(nghbr_counts);
			free_device(cumu_nghbr_counts);
		}
	}

} AssembledSolution;

}
}
}