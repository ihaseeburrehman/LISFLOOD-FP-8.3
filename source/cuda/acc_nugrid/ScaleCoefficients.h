#pragma once

#include "cuda_utils.cuh"
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct ScaleCoefficients
{
	NUMERIC_TYPE* h;
	NUMERIC_TYPE* v;
//	NUMERIC_TYPE* qx;
//	NUMERIC_TYPE* qy;
//	NUMERIC_TYPE* q3;
//	NUMERIC_TYPE* q4;
	NUMERIC_TYPE* z0;
	NUMERIC_TYPE* z1x;
	NUMERIC_TYPE* z1y;
//	NUMERIC_TYPE* zxy;
	NUMERIC_TYPE* n0;
	bool  is_copy = false;

	ScaleCoefficients(const int& num_all_elems, bool non_uniform_n)
	{
		size_t bytes = sizeof(NUMERIC_TYPE) * num_all_elems;
		
		h = (NUMERIC_TYPE*)malloc_device(bytes);
		v = (NUMERIC_TYPE*)malloc_device(bytes);
//		qx  = (NUMERIC_TYPE*)malloc_device(bytes);
//		qy  = (NUMERIC_TYPE*)malloc_device(bytes);
//		q3   = (NUMERIC_TYPE*)malloc_device(bytes);
//		q4 = (NUMERIC_TYPE*)malloc_device(bytes);
		z0 = (NUMERIC_TYPE*)malloc_device(bytes);
		z1x = (NUMERIC_TYPE*)malloc_device(bytes);
		z1y = (NUMERIC_TYPE*)malloc_device(bytes);
//		zxy = (NUMERIC_TYPE*)malloc_device(bytes);

		n0 = (non_uniform_n) ? (NUMERIC_TYPE*)malloc_device(bytes) : nullptr;

//		is_copy = false;
	}

	ScaleCoefficients(const ScaleCoefficients& original) { *this = original; is_copy = true; }

	~ScaleCoefficients()
	{
		if (!is_copy)
		{
			free_device(h);
			free_device(v);
//			free_device(qx);
//			free_device(qy);
//			free_device(q3);
//			free_device(q4);
			free_device(z0);
			free_device(z1x);
			free_device(z1y);
//			free_device(zxy);
			free_device(n0);
		}
	}

} ScaleCoefficients;

}
}
}