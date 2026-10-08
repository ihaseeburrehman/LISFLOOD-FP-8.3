#include "find_nonuniform_neighbours.cuh"

__global__ void lis::cuda::acc_nugrid::find_nonuniform_neighbours
(
	AssembledSolution    d_assem_sol, // compacted z order
	AssembledSolution    d_nghbr_assem_sol, // non-compacted row major
	Pars pars,
	Solver solver,
	NonUniformNeighbours d_non_uniform_nghbrs,
	Boundaries           boundaries
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_assem_sol.length) return;

	int level = d_assem_sol.levels[idx];
	
	index_1D act_idx = d_assem_sol.act_idcs[idx];

	MortonCode code = act_idx - get_lvl_idx(level);
	MortonCode code_sw = code << (2 * (solver.L - level));
	MortonCode code_ne = code << (2 * (solver.L - level)) ^ (0xFFFFFFFF >> (32 - 2 * (solver.L - level)));

	int x_sw = compact(code_sw);
	int y_sw = compact(code_sw >> 1);

	//if (x_sw == 0 && y_sw == 0) printf("%f\n", d_assem_sol.z0[idx]);

	int x_ne = compact(code_ne);
	int y_ne = compact(code_ne >> 1);

	int side_len = 1 << (solver.L - level);

	int mesh_dim = 1 << solver.L;

	int cumu_nghbr_count = d_assem_sol.cumu_nghbr_counts[idx];
	int count = 0;

	index_1D previous_nghbr = -5;

	// northern nghbrs using northeast point
	for (int i = 0; i < side_len; i++)
	{
		bool at_border = ( y_ne == (pars.ysz - 1) || ( y_ne == (mesh_dim - 1) ) );

		int x = x_ne - i;

		bool bc_bound = boundaries.north.bound(x);

		index_1D nghbr = (at_border) ? I_NORTH : (y_ne + 1) * mesh_dim + x;

		index_1D nghbr_check = (nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];

		bool already_in_list = (nghbr_check == previous_nghbr);

		previous_nghbr = nghbr_check;

		if (!already_in_list || (bc_bound && at_border))
		{
			d_non_uniform_nghbrs.nghbr_act_idcs[cumu_nghbr_count] = nghbr_check; //(nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];
			d_non_uniform_nghbrs.dx[cumu_nghbr_count] = (nghbr < 0) ? pars.dx * (1 << (solver.L - level)) : pars.dx * (1 << (solver.L - d_nghbr_assem_sol.levels[nghbr]));
			d_non_uniform_nghbrs.dirs[cumu_nghbr_count] = I_NORTH;
			cumu_nghbr_count++;
			count++;
		}
	}

	previous_nghbr = -5;

	// eastern nghbrs using northeast point
	for (int i = 0; i < side_len; i++)
	{
		bool at_border = ( x_ne == (pars.xsz - 1) || ( x_ne == (mesh_dim - 1) ) );

		int y = y_ne - i;

		bool bc_bound = boundaries.east.bound(y);

		index_1D nghbr = (at_border) ? I_EAST : y * mesh_dim + (x_ne + 1);

		index_1D nghbr_check = (nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];

		bool already_in_list = (nghbr_check == previous_nghbr);

		previous_nghbr = nghbr_check;

		//if (x_ne == 161 && y_ne == 0)
		//{
		//	printf
		//	(
		//		"x_ne: %d\n"
		//		"y_ne: %d\n"
		//		"y: %d\n"
		//		"ngbhr: %d\n"
		//		"nghbr_check: %d\n"
		//		"BC bound: %s\n"
		//		"In list: %s\n"
		//		"Border: %s\n"
		//		"Direction: %d\n"
		//		,
		//		x_ne, y_ne, y, nghbr, nghbr_check,
		//		bc_bound ? "YES" : "NO",
		//		already_in_list ? "YES" : "NO",
		//		at_border ? "YES" : "NO",
		//		I_EAST
		//	);
		//}

		if (!already_in_list  || (bc_bound && at_border) )
		{
			d_non_uniform_nghbrs.nghbr_act_idcs[cumu_nghbr_count] = nghbr_check; //(nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];
			d_non_uniform_nghbrs.dx[cumu_nghbr_count] = (nghbr < 0) ? pars.dx * (1 << (solver.L - level)) : pars.dx * (1 << (solver.L - d_nghbr_assem_sol.levels[nghbr]));
			d_non_uniform_nghbrs.dirs[cumu_nghbr_count] = I_EAST;
			cumu_nghbr_count++;
			count++;
		}
	}

	previous_nghbr = -5;

	// southern nghbrs using southwest point
	for (int i = 0; i < side_len; i++)
	{
		bool at_border = (y_sw == 0);

		int x = x_sw + i;

		bool bc_bound = boundaries.south.bound(x);

		index_1D nghbr = (at_border) ? I_SOUTH : (y_sw - 1) * mesh_dim + x;

		index_1D nghbr_check = (nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];

		bool already_in_list = (nghbr_check == previous_nghbr);

		previous_nghbr = nghbr_check;

		if (!already_in_list  || (bc_bound && at_border))
		{
			d_non_uniform_nghbrs.nghbr_act_idcs[cumu_nghbr_count] = nghbr_check; //(nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];
			d_non_uniform_nghbrs.dx[cumu_nghbr_count] = (nghbr < 0) ? pars.dx * (1 << (solver.L - level)) : pars.dx * (1 << (solver.L - d_nghbr_assem_sol.levels[nghbr]));
			d_non_uniform_nghbrs.dirs[cumu_nghbr_count] = I_SOUTH;
			cumu_nghbr_count++;
			count++;
		}
	}

	previous_nghbr = -5;

	// western nghbrs using southwest point
	for (int i = 0; i < side_len; i++)
	{
		bool at_border = (x_sw == 0);

		int y = y_sw + i;

		bool bc_bound = boundaries.west.bound(y);

		index_1D nghbr = (at_border) ? I_WEST : y * mesh_dim + (x_sw - 1);

		index_1D nghbr_check = (nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];

		bool already_in_list = (nghbr_check == previous_nghbr);

		previous_nghbr = nghbr_check;

		if (!already_in_list  || (bc_bound && at_border))
		{
			d_non_uniform_nghbrs.nghbr_act_idcs[cumu_nghbr_count] = nghbr_check; //(nghbr < 0) ? nghbr : d_nghbr_assem_sol.act_idcs[nghbr];
			d_non_uniform_nghbrs.dx[cumu_nghbr_count] = (nghbr < 0) ? pars.dx * (1 << (solver.L - level)) : pars.dx * (1 << (solver.L - d_nghbr_assem_sol.levels[nghbr]));
			d_non_uniform_nghbrs.dirs[cumu_nghbr_count] = I_WEST;
			cumu_nghbr_count++;
			count++;
		}
	}

	int owner = idx;

	int nghbr_count = d_assem_sol.nghbr_counts[idx];

	idx = d_assem_sol.cumu_nghbr_counts[idx];

	for (int i = 0; i < nghbr_count; i++)
	{
		d_non_uniform_nghbrs.owner_elem_idcs[idx + i] = owner;
		d_non_uniform_nghbrs.owner_counts[idx + i] = count;
	}
}