#include "load_interface_q_vol.cuh"
#include <algorithm>

__global__ void lis::cuda::acc_nugrid::load_interface_q_vol
(
	NonUniformNeighbours d_non_uniform_nghbrs,
	NonUniformInterfaces d_non_uniform_itrfaces
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_non_uniform_itrfaces.length) return;

	d_non_uniform_itrfaces.q_vol[idx] = d_non_uniform_nghbrs.q[d_non_uniform_itrfaces.load_idcs[idx]]; // from previous time step
}