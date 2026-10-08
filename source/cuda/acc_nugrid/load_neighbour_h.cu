#include "load_neighbour_h.cuh"
#include <algorithm>

__global__ void lis::cuda::acc_nugrid::load_neighbour_h
(
	AssembledSolution d_assem_sol,
	NonUniformNeighbours d_non_uniform_nghbrs,
	bool non_uniform_n
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_non_uniform_nghbrs.length) return;

	d_non_uniform_nghbrs.h_owner[idx] = d_assem_sol.h[d_non_uniform_nghbrs.owner_elem_idcs[idx]];
	d_non_uniform_nghbrs.h_nghbr[idx] = d_assem_sol.h[d_non_uniform_nghbrs.nghbr_elem_idcs[idx]];
	d_non_uniform_nghbrs.z_owner[idx] = d_assem_sol.z0[d_non_uniform_nghbrs.owner_elem_idcs[idx]];
	d_non_uniform_nghbrs.z_nghbr[idx] = d_assem_sol.z0[d_non_uniform_nghbrs.nghbr_elem_idcs[idx]];

	if (non_uniform_n) {
		d_non_uniform_nghbrs.n_owner[idx] = d_assem_sol.n0[d_non_uniform_nghbrs.owner_elem_idcs[idx]];
		d_non_uniform_nghbrs.n_nghbr[idx] = d_assem_sol.n0[d_non_uniform_nghbrs.nghbr_elem_idcs[idx]];
	}
}