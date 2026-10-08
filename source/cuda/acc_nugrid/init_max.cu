#include "init_max.cuh"
#include <algorithm>

__global__ void lis::cuda::acc_nugrid::init_max
(
	AssembledSolution    d_assem_sol,
	NUMERIC_TYPE*        maxH,
	NUMERIC_TYPE*        totalHtm,
	NUMERIC_TYPE*        maxHtm,
	NUMERIC_TYPE*        initHtm
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_assem_sol.length) return;

	maxHtm[idx] = NULLVAL; 
	initHtm[idx] = NULLVAL; 
	maxH[idx] = C(0.0);
	totalHtm[idx] = C(0.0);

}
