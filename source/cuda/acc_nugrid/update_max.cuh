#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "AssembledSolution.h"
#include "Directions.h"
#include "MortonCode.h"
#include "get_lvl_idx.cuh"
#include "get_level.cuh"
#include "NonUniformNeighbours.h"
#include "NonUniformInterfaces.h"
#include "GhostCellTypes.h"
#include "compact.cuh"
#include "Boundaries.h"
#include "PointSources.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void update_max
(
	AssembledSolution    d_assem_sol,
	Solver solver,
	NUMERIC_TYPE* maxH,
	NUMERIC_TYPE* totalHtm,
	NUMERIC_TYPE* maxHtm,
	NUMERIC_TYPE* initHtm
);

}
}
}