#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "ScaleCoefficients.h"
#include "AssembledSolution.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void reinsert_assem_sol
(
	AssembledSolution d_assem_sol, 
	ScaleCoefficients d_scale_coeffs
);

}
}
}