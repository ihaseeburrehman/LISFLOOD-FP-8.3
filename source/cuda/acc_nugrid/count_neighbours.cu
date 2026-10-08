#include "count_neighbours.cuh"

__global__ void lis::cuda::acc_nugrid::count_neighbours
(
    AssembledSolution d_assem_sol,
	AssembledSolution d_nghbr_assem_sol,
	Pars pars,
	Solver solver,
    Boundaries boundaries
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_assem_sol.length) return;

	int level = d_assem_sol.levels[idx];

	MortonCode code = d_assem_sol.act_idcs[idx] - get_lvl_idx(level); // MKS: code in level
	MortonCode code_sw = code << (2 * (solver.L - level));
	MortonCode code_ne = code << (2 * (solver.L - level)) ^ (0xFFFFFFFF >> (32 - 2 * (solver.L - level)));

	int x_sw = compact(code_sw);
	int y_sw = compact(code_sw >> 1);

	int x_ne = compact(code_ne);
	int y_ne = compact(code_ne >> 1);

	int mesh_dim = 1 << solver.L;

	int side_len = 1 << (solver.L - level);

//	index_1D all_nghbrs[4 * side_len]; // in order N, E, S, W

	index_1D previous_nghbr = -5;

	int nghbr_count = 0;

	// northern nghbrs using northeast point
	for (int i = 0; i < side_len; i++)
	{
		bool at_border = (y_ne == (pars.ysz - 1));

		int x = x_ne - i;

		bool bc_bound = boundaries.north.bound(x);

		index_1D nghbr = (at_border) ? I_NORTH : (y_ne + 1) * mesh_dim + x; // at_border-> nghbr = -1

		index_1D nghbr_ckeck = (nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];

		bool already_in_list = (nghbr_ckeck == previous_nghbr);

		previous_nghbr = nghbr_ckeck;

		if (!already_in_list || (bc_bound && at_border)) nghbr_count++;
	}

	previous_nghbr = -5;

	// eastern nghbrs using northeast point
	for (int i = 0; i < side_len; i++)
	{
		bool at_border = (x_ne == (pars.xsz - 1));

		int y = y_ne - i;

		bool bc_bound = boundaries.east.bound(y);

		index_1D nghbr = (at_border) ? I_EAST : y * mesh_dim + (x_ne + 1); // at_border-> nghbr = -2

		index_1D nghbr_ckeck = (nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];

		bool already_in_list = (nghbr_ckeck == previous_nghbr);

		previous_nghbr = nghbr_ckeck;

		if (!already_in_list || (bc_bound && at_border)) nghbr_count++;

	}

	previous_nghbr = -5;

	// southern nghbrs using southwest point
	for (int i = 0; i < side_len; i++)
	{
		bool at_border = (y_sw == 0);

		int x = x_sw + i;

		bool bc_bound = boundaries.south.bound(x);

		index_1D nghbr = (at_border) ? I_SOUTH : (y_sw - 1) * mesh_dim + x; // at_border-> nghbr = -3

		index_1D nghbr_ckeck = (nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];

		bool already_in_list = (nghbr_ckeck == previous_nghbr);

		previous_nghbr = nghbr_ckeck;

		if (!already_in_list || (bc_bound && at_border)) nghbr_count++;

	}

	previous_nghbr = -5;

	// western nghbrs using southwest point
	for (int i = 0; i < side_len; i++)
	{
		bool at_border = (x_sw == 0);

		int y = y_sw + i;

		bool bc_bound = boundaries.west.bound(y);

		index_1D nghbr = (at_border) ? I_WEST : y * mesh_dim + (x_sw - 1); // at_border-> nghbr = -4

		index_1D nghbr_ckeck = (nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];

		bool already_in_list = (nghbr_ckeck == previous_nghbr);

		previous_nghbr = nghbr_ckeck;

		if (!already_in_list || (bc_bound && at_border)) nghbr_count++;

	}

	d_assem_sol.nghbr_counts[idx] = nghbr_count;

}