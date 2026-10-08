#include "init_q.cuh"
#include <algorithm>

__global__ void lis::cuda::acc_nugrid::init_q
(
	NonUniformNeighbours    d_non_uniform_nghbrs
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_non_uniform_nghbrs.length) return;

	d_non_uniform_nghbrs.q[idx] = C(0.0);
	d_non_uniform_nghbrs.v[idx] = C(0.0);

}