#pragma once

#include "cuda_runtime.h"
#include "ScaleCoefficients.h"
#include "ParentScaleCoefficient.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__device__ __forceinline__ ParentScaleCoefficient load_parent_scale_coefficients
(
	ScaleCoefficients& d_scale_coeffs,
	int&               g_idx
)
{
	ParentScaleCoefficient parent_scale_coeffs =
	{
		d_scale_coeffs.h[g_idx],
		//d_scale_coeffs.qx[g_idx],
		//d_scale_coeffs.qy[g_idx],
		//d_scale_coeffs.z[g_idx]

		d_scale_coeffs.z0[g_idx],
		d_scale_coeffs.z1x[g_idx],
		d_scale_coeffs.z1y[g_idx],
//		d_scale_coeffs.zxy[g_idx]

	};

	return parent_scale_coeffs;
}

}
}
}