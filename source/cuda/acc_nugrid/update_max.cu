#include "update_max.cuh"
#include <algorithm>

__global__ void lis::cuda::acc_nugrid::update_max
(
	AssembledSolution    d_assem_sol,
	Solver solver,
	NUMERIC_TYPE*        maxH,
	NUMERIC_TYPE*        totalHtm,
	NUMERIC_TYPE*        maxHtm,
	NUMERIC_TYPE*        initHtm
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_assem_sol.length) return;

	if (d_assem_sol.h[idx] > solver.DepthThresh)
	{
		if (initHtm[idx] == NULLVAL)
			initHtm[idx] = solver.t / C(3600.0);

		totalHtm[idx] += solver.Tstep / C(3600.0);
		
		if (d_assem_sol.h[idx] > maxH[idx])
		{
			maxH[idx] = d_assem_sol.h[idx];
			maxHtm[idx] = solver.t / C(3600.0);
		}
	}
}
