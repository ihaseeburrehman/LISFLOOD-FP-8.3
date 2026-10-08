#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "AssembledSolution.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void 
calculate_dt
(
    AssembledSolution    d_assem_sol,
    Pars pars,
    Solver solver,
    NUMERIC_TYPE*                d_dt_CFL
);

}
}
}