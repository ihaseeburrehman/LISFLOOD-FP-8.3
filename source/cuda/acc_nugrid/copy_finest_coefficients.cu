#include "copy_finest_coefficients.cuh"

__global__ void lis::cuda::acc_nugrid::copy_finest_coefficients
(
	AssembledSolution d_assembled_solution,
	ScaleCoefficients d_scale_coeffs,
	index_1D          finest_lvl_idx,
	bool              non_uniform_n
)
{
	index_1D idx = blockDim.x * blockIdx.x + threadIdx.x;

	if ( idx >= d_assembled_solution.length ) return;
	
	index_1D g_idx = finest_lvl_idx + idx;

	d_scale_coeffs.h[g_idx] = d_assembled_solution.h[idx];
	//d_scale_coeffs.q1[g_idx] = d_assembled_solution.q1[idx];
	//d_scale_coeffs.q2[g_idx] = d_assembled_solution.q2[idx];
	//d_scale_coeffs.q3[g_idx] = d_assembled_solution.q3[idx];
	//d_scale_coeffs.q4[g_idx] = d_assembled_solution.q4[idx];

	d_scale_coeffs.z0[g_idx]  = d_assembled_solution.z0[idx];
	d_scale_coeffs.z1x[g_idx] = d_assembled_solution.z1x[idx];
	d_scale_coeffs.z1y[g_idx] = d_assembled_solution.z1y[idx];

	if (non_uniform_n) {
		d_scale_coeffs.n0[g_idx] = d_assembled_solution.n0[idx];
	}
}