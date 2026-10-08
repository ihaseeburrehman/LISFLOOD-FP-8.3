#include "count_interfaces_per_neighbours.cuh"

__global__ void lis::cuda::acc_nugrid::count_interfaces_per_neighbours
(
	NonUniformNeighbours d_non_uniform_nghbrs,
	AssembledSolution d_assem_sol
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_non_uniform_nghbrs.length) return;

	int interface_count = 0;
	
	int nghbr_dir = d_non_uniform_nghbrs.dirs[idx];

	// FINDING OWNER INTERFACES //

	int owner_elem_idx = d_non_uniform_nghbrs.owner_elem_idcs[idx];

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

		if (vertical || horizontal) interface_count++;



	}

	// ------------------------ //

	index_1D nghbr_act_idx = d_non_uniform_nghbrs.nghbr_act_idcs[idx];



	// if at ghost neighbour, return neighbours of the owner plus 2 as
	// ghost cell is at equal resolution with two interfaces
	if (nghbr_act_idx < 0)
	{
	//	d_non_uniform_nghbrs.itrface_counts[idx] = interface_count + 2;
		return;
	}

	// FINDING NEIGHBOUR INTERFACES //

	int nghbr_elem_idx = 0;

	for (int i = 0; i < d_assem_sol.length; i++)
	{
		nghbr_elem_idx = i;

		index_1D soln_act_idx = d_assem_sol.act_idcs[i];

		if (nghbr_act_idx == soln_act_idx) break;
	}

	// saving neighbour element index with which to access assem sol 
	d_non_uniform_nghbrs.nghbr_elem_idcs[idx] = nghbr_elem_idx;

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

		if (vertical || horizontal) interface_count++;
	}

	// ---------------------------- //

	d_non_uniform_nghbrs.itrface_counts[idx] = interface_count;

	// SETTING FLAGS FOR h0, z0, h1, z1 //

	if (nghbr_dir == I_EAST || nghbr_dir ==  I_SOUTH)
	{
		d_non_uniform_nghbrs.owner_flags[idx] = 1;
		d_non_uniform_nghbrs.nghbr_flags[idx] = -1;
	}
	else if (nghbr_dir == I_WEST || nghbr_dir == I_NORTH)
	{
		d_non_uniform_nghbrs.owner_flags[idx] = -1;
		d_non_uniform_nghbrs.nghbr_flags[idx] = 1;
	}

	// -------------------------------- //
}