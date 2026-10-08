#pragma once

#include "cuda_utils.cuh"
#include "encode_and_thresh_topo.cuh"
#include "get_num_blocks.h"
#include "index_1D.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ void preflag_topo
(
	ScaleCoefficients& d_scale_coeffs,
	Details&           d_details,
	bool*              d_sig_details,
	Maxes&             maxes, 
	NUMERIC_TYPE&      eps,
	int&               lev,
	bool               non_uniform_n,
	int                startfile
);

}
}
}