#pragma once

#include "../../rain/rain.h"
#include "BLOCK_VAR_MACROS.cuh"
#include "AssembledSolution.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "get_lvl_idx.cuh"
#include "MortonCode.h"
#include "compact.cuh"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__device__
NUMERIC_TYPE rate_at_cell
(
	NUMERIC_TYPE* rain_data,
	int i,
	int j,
	NUMERIC_TYPE dx,
	Pars pars,
	int xsz_rain,
	NUMERIC_TYPE tly_rain,
	NUMERIC_TYPE dx_rain,
	NUMERIC_TYPE dy_rain
);
	
//namespace DynamicRain
//{
__global__
void update_rain
(
	AssembledSolution    d_assem_sol,
	NUMERIC_TYPE* rain_data,
	int xsz_rain,
	NUMERIC_TYPE tly_rain,
	NUMERIC_TYPE dx_rain,
	NUMERIC_TYPE dy_rain,
	Pars pars,
	Solver solver
);

//}

__global__
void drain_nodata_water
(
	AssembledSolution    d_assem_sol,
	Pars pars,
	Solver solver
);

}
}
}