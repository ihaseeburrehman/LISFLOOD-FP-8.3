#pragma once

#include <stdio.h>
#include <string.h>
#include "BLOCK_VAR_MACROS.cuh"
#include "cuda_utils.cuh"
#include "AssembledSolution.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ void write_raster_maps
(
	const char*                 respath,
	const AssembledSolution&    d_assem_sol,
	const int&                  mesh_dim,
	const Pars& pars,
	const int& call_gzip,
	const int& precision = DEFAULT_PRECISION
);

}
}
}