#pragma once

#include "cuda_runtime.h"
#include "Detail.h"
#include "Details.h"
#include "index_1D.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__device__ __forceinline__ Detail load_details
(
	const Details&  d_details,
	const index_1D& g_idx
)
{
	Detail detail;

	//detail.eta.alpha = d_details.eta.alpha[g_idx];
	//detail.eta.beta = d_details.eta.beta[g_idx];
	//detail.eta.gamma = d_details.eta.gamma[g_idx];

	//detail.qx.alpha = d_details.qx.alpha[g_idx];
	//detail.qx.beta = d_details.qx.beta[g_idx];
	//detail.qx.gamma = d_details.qx.gamma[g_idx];

	//detail.qy.alpha = d_details.qy.alpha[g_idx];
	//detail.qy.beta = d_details.qy.beta[g_idx];
	//detail.qy.gamma = d_details.qy.gamma[g_idx];

	detail.z.alpha_11 = d_details.z0.alpha[g_idx];
	detail.z.beta_11 = d_details.z0.beta[g_idx];
	detail.z.gamma_11 = d_details.z0.gamma[g_idx];

	detail.z.alpha_12 = d_details.z1y.alpha[g_idx];
	detail.z.beta_12 = d_details.z1y.beta[g_idx];
	detail.z.gamma_12 = d_details.z1y.gamma[g_idx];

	detail.z.alpha_21 = d_details.z1x.alpha[g_idx];
	detail.z.beta_21 = d_details.z1x.beta[g_idx];
	detail.z.gamma_21 = d_details.z1x.gamma[g_idx];

//	detail.z.alpha_22 = d_details.zxy.alpha[g_idx];
//	detail.z.beta_22 = d_details.zxy.beta[g_idx];
//	detail.z.gamma_22 = d_details.zxy.gamma[g_idx];

	return detail;
}

}
}
}