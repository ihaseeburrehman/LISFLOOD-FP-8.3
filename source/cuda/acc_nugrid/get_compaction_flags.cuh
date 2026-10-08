#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cub/block/block_load.cuh"
#include "BLOCK_VAR_MACROS.cuh"
#include "index_1D.h"
#include "CompactionFlags.h"
#include "AssembledSolution.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void get_compaction_flags
(
	AssembledSolution d_solution,
	CompactionFlags   d_compaction_flags,
	int               num_finest_elems
);

}
}
}