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

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ void get_modal_values
(
	AssembledSolution& d_assembled_solution,
	const int&         mesh_dim,
	const int&         interface_dim,
	const Fnames&      filenames,
	const States&      states,
	const NUMERIC_TYPE& no_data
);

}
}
}