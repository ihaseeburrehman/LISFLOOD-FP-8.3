#pragma once

#include "cuda_utils.cuh"
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct SubDetails
{
	NUMERIC_TYPE* alpha;
	NUMERIC_TYPE* beta;
	NUMERIC_TYPE* gamma;
	bool  is_copy = false;

	SubDetails() = default; //new

	SubDetails(const int& num_details)
	{
		size_t bytes = sizeof(NUMERIC_TYPE) * num_details;

		alpha = (NUMERIC_TYPE*)malloc_device(bytes);
		beta  = (NUMERIC_TYPE*)malloc_device(bytes);
		gamma = (NUMERIC_TYPE*)malloc_device(bytes);

//		is_copy = false;
	}

	SubDetails(const SubDetails& original) { *this = original; is_copy = true; }

	~SubDetails()
	{
		if (!is_copy)
		{
			free_device(alpha);
			free_device(beta);
			free_device(gamma);
		}
	}

} SubDetails;

}
}
}