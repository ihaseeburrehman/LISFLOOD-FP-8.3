#pragma once

#ifndef __CUDACC__
    #define __CUDACC__
#endif

#include "cuda_runtime.h"
#include "cub/cub.cuh"
#include "CHECK_CUDA_ERROR.cuh"
#include "cuda_utils.cuh"
#include "AssembledSolution.h"
#include "MortonCode.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ void sort_finest_scale_coefficients_z_order
(
	MortonCode*        d_morton_codes,
	MortonCode*        d_sorted_morton_codes,
	AssembledSolution& d_assembled_solution,
	AssembledSolution& d_buffer_assembled_solution,
	MortonCode*        d_indices,
	MortonCode*        d_reverse_z_order,
	bool               non_uniform_n
);

}
}
}