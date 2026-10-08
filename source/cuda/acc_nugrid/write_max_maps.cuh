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

__host__ void write_max_maps
(
	const char*                 respath,
	const int&                  mesh_dim,
	const Pars& pars,
	NUMERIC_TYPE* maxH,
	NUMERIC_TYPE* totalHtm,
	NUMERIC_TYPE* maxHtm,
	NUMERIC_TYPE* initHtm,
	const int& precision = DEFAULT_PRECISION
);

}
}
}