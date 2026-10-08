#pragma once

#include "cub/device/device_select.cuh"

#include "CHECK_CUDA_ERROR.cuh"
#include "cuda_utils.cuh"
#include "AssembledSolution.h"
#include "CompactionFlags.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

void compaction
(
	AssembledSolution& d_sol, 
	AssembledSolution& d_buf_sol, 
	CompactionFlags&   d_compaction_flags,
	int                num_finest_elems,
	bool               non_uniform_n,
	int                startfile
);

}
}
}