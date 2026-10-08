#pragma once

#include "cuda_runtime.h"
#include "CHECK_CUDA_ERROR.cuh"
#include "BLOCK_VAR_MACROS.cuh"
#include "cuda_utils.cuh"
#include "reg_new.cuh"
#include "get_num_blocks.h"
#include "index_1D.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ void get_reg_tree
(
	bool*            d_sig_details,
	int           lev
);

}
}
}