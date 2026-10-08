#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "Directions.h"
#include "MortonCode.h"
#include "get_lvl_idx.cuh"
#include "get_level.cuh"
#include "NonUniformNeighbours.h"
#include "GhostCellTypes.h"
#include "compact.cuh"
#include "Boundaries.h"
#include "NonUniformInterfaces.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void load_interface_q_vol
(
	NonUniformNeighbours d_non_uniform_nghbrs,
	NonUniformInterfaces d_non_uniform_itrfaces
);

}
}
}