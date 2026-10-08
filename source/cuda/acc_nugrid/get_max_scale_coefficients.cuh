#pragma once

#include "cuda_runtime.h"
#include <cmath>
#include "CHECK_CUDA_ERROR.cuh"
#include "BLOCK_VAR_MACROS.cuh"
#include "cuda_utils.cuh"
#include "AssembledSolution.h"
#include "Maxes.h"
#include "get_num_blocks.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

struct AbsMax
{
	template <typename T>
	__device__ __forceinline__
	T operator()(const T& a, const T& b) const { return ( FABS(a) > FABS(b) ) ? FABS(a) : FABS(b); }
};

__host__ Maxes get_max_scale_coefficients
(
	AssembledSolution& d_assem_sol
);

}
}
}