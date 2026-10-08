#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "../../lisflood.h"
#include "DetailChildren.h"
#include "BLOCK_VAR_MACROS.cuh"
#include "index_1D.h"
#include "get_child_details.cuh"
#include "get_lvl_idx.cuh"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

template <bool SINGLE_BLOCK>
__global__ void reg_new
(
	bool*    d_sig_details,
	int      level,
	index_1D prev_lvl_idx,
	index_1D curr_lvl_idx,
	index_1D next_lvl_idx,
	int      num_threads
);

}
}
}