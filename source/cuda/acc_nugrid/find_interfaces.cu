#include "find_interfaces.cuh"
#include <algorithm>

__global__ void lis::cuda::acc_nugrid::find_interfaces
(
	NonUniformNeighbours d_non_uniform_nghbrs,
	AssembledSolution d_assem_sol,
	NonUniformInterfaces d_non_uniform_itrfaces
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_non_uniform_nghbrs.length) return;

	int cumu_interface_count = d_non_uniform_nghbrs.cumu_itrface_counts[idx];

	int nghbr_dir = d_non_uniform_nghbrs.dirs[idx]; // MKS

	NUMERIC_TYPE dx_sum = C(0.0);
	NUMERIC_TYPE dx     = C(0.0);

	int count = 0;

	// FINDING OWNER INTERFACES //

	int owner_elem_idx = d_non_uniform_nghbrs.owner_elem_idcs[idx];

	//
	//int level = d_assem_sol.levels[d_non_uniform_nghbrs.owner_elem_idcs[idx]];
	//MortonCode code = d_assem_sol.act_idcs[d_non_uniform_nghbrs.owner_elem_idcs[idx]] - get_lvl_idx(level);
	//int x = compact(code);
	//int y = compact(code >> 1);
	//

	for (int i = 0; i < d_assem_sol.nghbr_counts[owner_elem_idx]; i++)
	{
		int owner_itrfaces_idx = d_assem_sol.cumu_nghbr_counts[owner_elem_idx] + i;

		int interface_dir = d_non_uniform_nghbrs.dirs[owner_itrfaces_idx];

		bool vertical =
			(
				(
					nghbr_dir == I_EAST
					||
					nghbr_dir == I_WEST
					)
				&&
				(
					interface_dir == I_NORTH
					||
					interface_dir == I_SOUTH
					)
				);

		bool horizontal =
			(
				(
					nghbr_dir == I_NORTH
					||
					nghbr_dir == I_SOUTH
					)
				&&
				(
					interface_dir == I_EAST
					||
					interface_dir == I_WEST
					)
				);

		// ------------------------ //

		index_1D nghbr_act_idx = d_non_uniform_nghbrs.nghbr_act_idcs[idx];

		// TODO: TREATMENT OF GHOST INTERFACES //
		if (nghbr_act_idx < 0)
		{
			//		dxs[count] = -9999;
			d_non_uniform_nghbrs.nghbr_elem_idcs[idx] = 0;
			return;
		}
		// ----------------------------------- //


		if (vertical || horizontal)
		{
			dx = d_non_uniform_nghbrs.dx[owner_itrfaces_idx];

			dx_sum += dx;

			d_non_uniform_itrfaces.dx[cumu_interface_count] = dx;
			d_non_uniform_itrfaces.load_idcs[cumu_interface_count] = owner_itrfaces_idx;
			cumu_interface_count++;
			count++;
		}
	}

	// FINDING NEIGHBOUR INTERFACES //
	int nghbr_elem_idx = d_non_uniform_nghbrs.nghbr_elem_idcs[idx];

	for (int i = 0; i < d_assem_sol.nghbr_counts[nghbr_elem_idx]; i++)
	{
		int nghbr_itrfaces_idx = d_assem_sol.cumu_nghbr_counts[nghbr_elem_idx] + i;

		int interface_dir = d_non_uniform_nghbrs.dirs[nghbr_itrfaces_idx];

		bool vertical =
			(
				(
					nghbr_dir == I_EAST
					||
					nghbr_dir == I_WEST
					)
				&&
				(
					interface_dir == I_NORTH
					||
					interface_dir == I_SOUTH
					)
				);

		bool horizontal =
			(
				(
					nghbr_dir == I_NORTH
					||
					nghbr_dir == I_SOUTH
					)
				&&
				(
					interface_dir == I_EAST
					||
					interface_dir == I_WEST
					)
				);

		if (vertical || horizontal)
		{
			dx = d_non_uniform_nghbrs.dx[nghbr_itrfaces_idx];

			dx_sum += dx;

			d_non_uniform_itrfaces.dx[cumu_interface_count] = dx;
			d_non_uniform_itrfaces.load_idcs[cumu_interface_count] = nghbr_itrfaces_idx;
			cumu_interface_count++;
			count++;
		}
	}

	// ---------------------------- //

	d_non_uniform_nghbrs.dx_sum[idx] = dx_sum;

	int nghbr = idx;

	int interface_count = d_non_uniform_nghbrs.itrface_counts[idx];

	idx = d_non_uniform_nghbrs.cumu_itrface_counts[idx];

	for (int i = 0; i < interface_count; i++)
	{
		d_non_uniform_itrfaces.q_vol[idx + i] = C(0.0);
		d_non_uniform_itrfaces.nghbrs[idx + i] = nghbr;
		d_non_uniform_itrfaces.nghbr_counts[idx + i] = count;
	}
}