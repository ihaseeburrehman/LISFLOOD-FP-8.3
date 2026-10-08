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
#include "GhostCellTypes.h"
#include "compact.cuh"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void count_neighbours
(
    AssembledSolution d_assem_sol,
    AssembledSolution d_nghbr_assem_sol,
    Pars pars,
    Solver solver,
    Boundaries boundaries
);

}
}
}