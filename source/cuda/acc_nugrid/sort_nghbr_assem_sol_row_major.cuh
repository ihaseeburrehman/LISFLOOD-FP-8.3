#pragma once

#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BLOCK_VAR_MACROS.cuh"
#include "CHECK_CUDA_ERROR.cuh"
#include "cuda_utils.cuh"
#include "AssembledSolution.h"
#include "get_lvl_idx.cuh"
#include "get_num_blocks.h"
#include "read_raster_file.h"
#include "write_raster_maps.cuh"
#include "MortonCode.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ void sort_nghbr_assem_sol_row_major
(
	MortonCode*        d_reverse_z_order,
	MortonCode*        d_indices,
	AssembledSolution& d_buf_assem_sol,
	AssembledSolution& d_nghbr_assem_sol,
	int                array_length
);

}
}
}