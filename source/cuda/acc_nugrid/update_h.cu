#include "update_h.cuh"
#include <algorithm>

__global__ void lis::cuda::acc_nugrid::update_h
(
	AssembledSolution    d_assem_sol,
	NonUniformNeighbours d_non_uniform_nghbrs,
	NonUniformInterfaces d_non_uniform_itrfaces,
	Pars pars,
	Solver solver,
	Boundaries           boundaries,
	PointSources      point_sources
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_assem_sol.length) return;

	NUMERIC_TYPE flag;
	NUMERIC_TYPE q = C(0.0);

	int index = d_assem_sol.cumu_nghbr_counts[idx];





	//int level = d_assem_sol.levels[d_non_uniform_nghbrs.owner_elem_idcs[idx]];

	//MortonCode code = d_assem_sol.act_idcs[idx] - get_lvl_idx(level);
	//int x = compact(code);
	//int y = compact(code >> 1);





	// Check for drying elements
	for (int i = 0; i < d_assem_sol.nghbr_counts[idx]; i++)
	{
		if (d_non_uniform_nghbrs.dirs[index + i] == I_NORTH) flag = 1.0;
		if (d_non_uniform_nghbrs.dirs[index + i] == I_EAST) flag = -1.0;
		if (d_non_uniform_nghbrs.dirs[index + i] == I_SOUTH) flag = -1.0;
		if (d_non_uniform_nghbrs.dirs[index + i] == I_WEST) flag = 1.0;

	//	if (idx == 9840)
	//	{
	//		//printf("%d\n", g_idx);
	//		//printf("%d\n", get_lvl_idx(params.L));
	//		//printf("%d\n", point_sources.d_codes[i]);
	//		printf("idx %d\n", idx);
	////		printf("H %f\n", H);
	//		printf("qacc %f\n", q);
	//		printf("q %f\n", d_non_uniform_nghbrs.q[index + i]);
	//		printf("dir %d\n", d_non_uniform_nghbrs.dirs[index + i]);
	//		printf("i %d\n", i);
	//		printf("index %d\n", index);
	//		//	//printf("%d\n", previous_nghbr);
	//		//	//printf("%d\n", bc_bound);
	//		//	//		for (int i = 0; i < nghbr_count; i++)
	//		//	//		{
	//		//	////			printf("%d\n", all_nghbrs[i]);
	//		//	//		}
	//	}



		//if (x == 79 && y == 99)
		//{
		//	printf
		//	(
		//		"i: %d\n"
		//		"q: %15.12f\n\n"
		//		"flag: %f\n\n"
		//		,
		//		i,d_non_uniform_nghbrs.q[index + i], flag
		//	);
		//}

		q += (flag * d_non_uniform_nghbrs.q[index + i]);
	}

	NUMERIC_TYPE dV = solver.Tstep * q;
	NUMERIC_TYPE dx = pars.dx * (1 << (solver.L - d_assem_sol.levels[idx]));
	NUMERIC_TYPE H = d_assem_sol.h[idx];
	NUMERIC_TYPE cv = H * (dx * dx);

	//if (x == 79 && y == 99)
	//{
	//	printf
	//	(
	//		"dV: %15.12f\n"
	//		"H: %15.12f\n\n"
	//		"cv: %15.12f\n\n"
	//		,
	//		dV, H, cv
	//	);
	//}

	if (H > solver.DepthThresh)
	{
		if (cv + dV < C(0.0))
		{
			NUMERIC_TYPE WDweight = -cv / dV; // C(0.5) for improved drying stability
			for (int i = 0; i < d_assem_sol.nghbr_counts[idx]; i++)
			{
				if (d_non_uniform_nghbrs.dirs[index + i] == I_WEST && d_non_uniform_nghbrs.q[index + i] < C(0.0)) d_non_uniform_nghbrs.q[index + i] *= WDweight;
				if (d_non_uniform_nghbrs.dirs[index + i] == I_EAST && d_non_uniform_nghbrs.q[index + i] > C(0.0)) d_non_uniform_nghbrs.q[index + i] *= WDweight;
				if (d_non_uniform_nghbrs.dirs[index + i] == I_SOUTH && d_non_uniform_nghbrs.q[index + i] > C(0.0)) d_non_uniform_nghbrs.q[index + i] *= WDweight;
				if (d_non_uniform_nghbrs.dirs[index + i] == I_NORTH && d_non_uniform_nghbrs.q[index + i] < C(0.0)) d_non_uniform_nghbrs.q[index + i] *= WDweight;
			}
		}
	}

	//if (x == 79 && y == 99)
	//{
	//	printf
	//	(
	//		"q: %15.12f\n"
	//		"dV: %15.12f\n"
	//		"H: %15.12f\n\n"
	//		,
	//		q, dV, H
	//	);
	//}



	// End of Check for drying elements


	// point BC goes here



	q = C(0.0);
	NUMERIC_TYPE v_north = C(0.0);
	NUMERIC_TYPE v_east = C(0.0);
	NUMERIC_TYPE v_south = C(0.0);
	NUMERIC_TYPE v_west = C(0.0);

	int north_counter = 0;
	int east_counter = 0;
	int south_counter = 0;
	int west_counter = 0;

	for (int i = 0; i < d_assem_sol.nghbr_counts[idx]; i++)
	{
		if (d_non_uniform_nghbrs.dirs[index + i] == I_NORTH) {
			flag = C(1.0);
			v_north += (flag * d_non_uniform_nghbrs.v[index + i]);
			north_counter = north_counter + 1;
		}
		if (d_non_uniform_nghbrs.dirs[index + i] == I_EAST) {
			flag = -C(1.0);
			v_east += (flag * d_non_uniform_nghbrs.v[index + i]);
			east_counter = east_counter + 1;
		}
		if (d_non_uniform_nghbrs.dirs[index + i] == I_SOUTH) {
			flag = -C(1.0);
			v_south += (flag * d_non_uniform_nghbrs.v[index + i]);
			south_counter = south_counter + 1;
		}
		if (d_non_uniform_nghbrs.dirs[index + i] == I_WEST) {
			flag = C(1.0);
			v_west += (flag * d_non_uniform_nghbrs.v[index + i]);
			west_counter = west_counter + 1;
		}


		q += (flag * d_non_uniform_nghbrs.q[index + i]);

		//if (x == 79 && y == 99)
		//{
		//	printf
		//	(
		//		"i: %d\n"
		//		"q: %15.12f\n"
		//		"flag: %15.12f\n\n"

		//		,
		//		i, d_non_uniform_nghbrs.q[index + i], flag
		//	);
		//}

	}


	NUMERIC_TYPE v_north_sum = v_north/ north_counter;
	NUMERIC_TYPE v_east_sum = v_east / east_counter;
	NUMERIC_TYPE v_south_sum = v_south / south_counter;
	NUMERIC_TYPE v_west_sum = v_west / west_counter;

	d_assem_sol.v[idx]= sqrt(pow(getmax(fabs(v_west_sum), fabs(v_east_sum)), 2) + pow(getmax(fabs(v_north_sum), fabs(v_south_sum)), 2));

	dV = solver.Tstep * q;
	H += dV / (dx * dx);
	if (H < C(0.0)) H = C(0.0);

	//if (x == 79 && y == 99)
	//{
	//	printf
	//	(
	//		"q: %15.12f\n"
	//		"dV: %15.12f\n"
	//		"H: %15.12f\n\n"
	//		,
	//		q, dV, H
	//	);
	//}

	//if (y == 0 || y == 1)
	//{
	//	H = C(0.0);
	//}

//	if (point_sources.num_srcs)
//	{
//		for (int i = 0; i < point_sources.num_srcs; i++)
//		{
//			g_idx = /*get_lvl_idx(params.L) +*/ point_sources.d_codes[i];
//			src_type = point_sources.d_src_types[i];
//			if (d_assem_sol.act_idcs[idx] == g_idx) {
//				if (src_type == HFIX || src_type == HVAR)
//				{
//					d_assem_sol.h[idx] = point_sources.d_srcs[i] - d_assem_sol.z0[idx];
//				}
////				else if (src_type == QFIX || src_type == QVAR)
////				{
////					d_assem_sol.h[idx] += point_sources.q_source(dt, dx, i); // coverts q to h
////				}
//			}
//		}
//	}






	// point BC goes here
	index_1D g_idx;
	int src_type;

	if (point_sources.num_srcs)
	{
		for (int i = 0; i < point_sources.num_srcs; i++)
		{
			g_idx = get_lvl_idx(solver.L) + point_sources.d_codes[i];
			src_type = point_sources.d_src_types[i];

			//			if (i == 0)
			//			{
							//printf("%d\n", g_idx);
							//printf("%d\n", g_idx);
				//			printf("%d\n", d_assem_sol.act_idcs[idx]);
							//printf("%d\n", idx);
							//printf("%d\n", interface_count);
							////printf("%d\n", vertical);
							////printf("%d\n", horizontal);
							//printf("%d\n", nghbr_dir);
							////printf("%d\n", nghbr_act_idx);
							////	printf("%d\n", y_ne);
							//	//printf("%d\n", previous_nghbr);
							//	//printf("%d\n", bc_bound);
							//	//		for (int i = 0; i < nghbr_count; i++)
							//	//		{
							//	////			printf("%d\n", all_nghbrs[i]);
							//	//		}
			//			}


			if (d_assem_sol.act_idcs[idx] == g_idx) {



				if (src_type == HFIX || src_type == HVAR)
				{
					H = point_sources.d_srcs[i] - d_assem_sol.z0[idx];
					if (H < C(0.0)) H = C(0.0);
				}
				else if (src_type == QFIX || src_type == QVAR)
				{
					//d_assem_sol.h[idx] += point_sources.q_source(dt, dx, i); // coverts q to h
					H += point_sources.q_src(solver.Tstep, dx, i); // coverts q to h
					if (H < C(0.0)) H = C(0.0);



				}
			}
		}
	}





	d_assem_sol.h[idx] = H;

}