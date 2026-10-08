#include "get_compaction_flags.cuh"

__global__ void lis::cuda::acc_nugrid::get_compaction_flags
(
	AssembledSolution d_solution,
	CompactionFlags   d_compaction_flags,
	int                num_finest_elems
)
{
	index_1D idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= num_finest_elems) return;

	index_1D current = d_solution.act_idcs[idx];

	if ( idx < (num_finest_elems - 1) )
	{
		index_1D right = d_solution.act_idcs[idx + 1];

		d_compaction_flags.north_east[idx] = !(current == right);
	}
	else
	{
		d_compaction_flags.north_east[idx] = 1;
	}

	if (idx > 0)
	{
		index_1D left = d_solution.act_idcs[idx - 1];

		d_compaction_flags.south_west[idx] = !(current == left);
	}
	else
	{
		d_compaction_flags.south_west[idx] = 1;
	}
}