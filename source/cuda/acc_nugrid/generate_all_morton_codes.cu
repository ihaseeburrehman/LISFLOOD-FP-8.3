#include "generate_all_morton_codes.cuh"

__global__ void lis::cuda::acc_nugrid::generate_all_morton_codes
(
	MortonCode* d_morton_codes,
	int*        d_indices, 
	int         mesh_dim
)
{
	index_1D idx = blockDim.x * blockIdx.x + threadIdx.x;

	if ( idx >= (mesh_dim * mesh_dim) ) return;

	int x = idx % mesh_dim;
	int y = idx / mesh_dim;

	d_indices[idx] = idx;

	d_morton_codes[idx] = generate_morton_code(x, y);
}