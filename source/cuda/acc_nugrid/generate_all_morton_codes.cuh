#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_utils.cuh"
#include "generate_morton_code.cuh"
#include "MortonCode.h"
#include "index_1D.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void generate_all_morton_codes
(
	MortonCode* d_morton_codes,
	int*        d_indices,
	int         mesh_dim
);

}
}
}