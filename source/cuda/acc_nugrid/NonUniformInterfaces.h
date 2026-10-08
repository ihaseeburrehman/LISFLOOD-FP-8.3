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

typedef struct NonUniformInterfaces
{
	NUMERIC_TYPE* q_vol;
	NUMERIC_TYPE* dx;
	NUMERIC_TYPE* dx_sum;
	int*  nghbrs;
	int*  nghbr_counts;
	int*  load_idcs;
	int   length;
	bool  is_copy;

	NonUniformInterfaces(const int& num_interfaces)
	{
		size_t bytes_real = num_interfaces * sizeof(NUMERIC_TYPE);
		size_t bytes_int  = num_interfaces * sizeof(int);

		q_vol             = (NUMERIC_TYPE*)malloc_device(bytes_real);
		dx                = (NUMERIC_TYPE*)malloc_device(bytes_real);
		dx_sum            = (NUMERIC_TYPE*)malloc_device(bytes_real);
		nghbrs            = (int*)malloc_device(bytes_int);
		nghbr_counts      = (int*)malloc_device(bytes_int);
		load_idcs         = (int*)malloc_device(bytes_int);
		length            = num_interfaces;
		is_copy           = false;
	}

	NonUniformInterfaces(const NonUniformInterfaces& original) { *this = original; is_copy = true; }

	~NonUniformInterfaces()
	{
		if (!is_copy)
	    {
		    free_device(q_vol);
		    free_device(dx);
		    free_device(dx_sum);
		    free_device(nghbrs);
		    free_device(nghbr_counts);
		    free_device(load_idcs);
	    }
	}
} NonUniformInterfaces;

}
}
}