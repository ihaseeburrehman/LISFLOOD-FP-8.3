#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cub/block/block_scan.cuh"

#include "BLOCK_VAR_MACROS.cuh"
#include "ScaleCoefficients.h"
#include "Details.h"
#include "Maxes.h"
#include "ChildScaleCoefficients.h"
#include "ParentScaleCoefficient.h"
#include "Detail.h"
#include "get_lvl_idx.cuh"
//#include "load_child_scale_coefficients.cuh"
#include "encode_scale.cuh"
#include "encode_details.cuh"
#include "index_1D.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__global__ void encode_and_thresh_topo
(
	ScaleCoefficients d_scale_coeffs,
	Details           d_details,
	bool*             d_sig_details,
	Maxes             maxes,
	NUMERIC_TYPE      epsilon_local,
	int               level,
	bool              non_uniform_n,
	int               startfile
);

}
}
}