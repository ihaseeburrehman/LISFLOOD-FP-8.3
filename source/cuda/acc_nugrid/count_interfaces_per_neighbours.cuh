#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "AssembledSolution.h"
#include "Directions.h"
#include "MortonCode.h"
#include "get_lvl_idx.cuh"
#include "get_level.cuh"
#include "Boundaries.h"
#include "NonUniformNeighbours.h"
#include "GhostCellTypes.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void count_interfaces_per_neighbours
(
	NonUniformNeighbours d_non_uniform_nghbrs,
	AssembledSolution d_assem_sol
);

}
}
}