#pragma once

#include "../rain/rain.h"
#include "cuda_solver.cuh"

namespace lis
{
namespace cuda
{

__device__ 
NUMERIC_TYPE rate_at_cell
(
	NUMERIC_TYPE* rain_data,
	int i,
	int j
);

__device__ 
NUMERIC_TYPE rate_at_cell_ACC
(
	NUMERIC_TYPE* rain_data,
	int i,
	int j
);

namespace DynamicRain
{
	__global__
	__launch_bounds__(CUDA_BLOCK_SIZE)
	void
	update
	(
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* rain_data
	);
}

namespace DynamicRain
{
	__global__
	__launch_bounds__(CUDA_BLOCK_SIZE)
	void
	updateACC
	(
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* rain_data
	);
}

}
}
