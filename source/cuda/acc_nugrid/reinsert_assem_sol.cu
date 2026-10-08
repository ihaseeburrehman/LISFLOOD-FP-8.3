#include "reinsert_assem_sol.cuh"

__global__ void lis::cuda::acc_nugrid::reinsert_assem_sol
(
	AssembledSolution d_assem_sol,
	ScaleCoefficients d_scale_coeffs
)
{
	index_1D idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_assem_sol.length) return;

	index_1D active_idx = d_assem_sol.act_idcs[idx];

	d_scale_coeffs.h[active_idx] = d_assem_sol.h[idx]; // +d_assem_sol.z[idx];
	d_scale_coeffs.v[active_idx] = d_assem_sol.v[idx];
	//d_scale_coeffs.qx[active_idx]  = d_assem_sol.qx[idx];
	//d_scale_coeffs.qy[active_idx]  = d_assem_sol.qy[idx];
	d_scale_coeffs.z0[active_idx]   = d_assem_sol.z0[idx];
//	d_scale_coeffs.z1x[active_idx] = d_assem_sol.z1x[idx];
//	d_scale_coeffs.z1y[active_idx] = d_assem_sol.z1y[idx];
//	d_scale_coeffs.zxy[active_idx] = d_assem_sol.zxy[idx];
}