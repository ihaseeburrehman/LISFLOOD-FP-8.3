#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "AssembledSolution.h"
#include "ScaleCoefficients.h"
#include "index_1D.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void copy_finest_coefficients
(
	AssembledSolution d_assembled_solution,
	ScaleCoefficients d_scale_coeffs,
	index_1D          finest_lvl_idx,
	bool              non_uniform_n
);

}
}
}