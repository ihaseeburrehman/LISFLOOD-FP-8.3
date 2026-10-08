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

__global__ void update_h
(
	AssembledSolution    d_assem_sol,
	NonUniformNeighbours d_non_uniform_nghbrs,
	NonUniformInterfaces d_non_uniform_itrfaces,
	Pars pars,
	Solver solver,
	Boundaries           boundaries,
	PointSources      point_sources
);

}
}
}