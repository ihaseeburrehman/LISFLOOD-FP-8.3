#include "compute_q.cuh"
#include <algorithm>

__global__ void lis::cuda::acc_nugrid::compute_q
(
	AssembledSolution    d_assem_sol,
	NonUniformNeighbours d_non_uniform_nghbrs,
	NonUniformInterfaces d_non_uniform_itrfaces,
	Pars pars,
	Solver solver,
	Boundaries           boundaries,
	bool non_uniform_n
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_non_uniform_nghbrs.length) return;

	index_1D nghbr_act_idx = d_non_uniform_nghbrs.nghbr_act_idcs[idx];

	NUMERIC_TYPE h_nghbr, h_owner, z_nghbr, z_owner, y_nghbr, y_owner; // , dx_nghbr, dx_owner;
	NUMERIC_TYPE qvect, q_opp, q_opp_avg, q_old;
	NUMERIC_TYPE dh, Sf, hflow, MaxHflow, theta, length, Manning, g;
	int sign, dir, edge;
	int boundary_type;
	NUMERIC_TYPE boundary_value = C(0.0);
	NUMERIC_TYPE h0, h1, z1;
	
	theta = C(1.0);
	length = d_non_uniform_nghbrs.dx[idx];

	int level = d_assem_sol.levels[d_non_uniform_nghbrs.owner_elem_idcs[idx]];

	MortonCode code = d_assem_sol.act_idcs[d_non_uniform_nghbrs.owner_elem_idcs[idx]] - get_lvl_idx(level);
	int x = compact(code);
	int y = compact(code >> 1);
	
//	bool finest = (level == pars.L);
//	bool bound = false;
//	bool border = false;
//	bool finebdy = false;

		
	g = solver.g;

	if (nghbr_act_idx == -1 ) {
		sign = -1;
		dir = 1;
		edge = 1;
		boundary_type = boundaries.north.bdytype;
		if (boundaries.north.bound(x) && level == solver.L) {
			boundary_value = boundaries.north.inlet;
		}
		else {
			boundary_type = CLOSED;
		}
	}
	else if (nghbr_act_idx == -2 ) {
		sign = 1;
		dir = 2;
		edge = 2;
		boundary_type = boundaries.east.bdytype;
		if (boundaries.east.bound(y) && level == solver.L) {
			boundary_value = boundaries.east.inlet;
		}
		else {
			boundary_type = CLOSED;
		}
	}
	else if (nghbr_act_idx == -3 ) {


		sign = 1;
		dir = 1;
		edge = 3;
		boundary_type = boundaries.south.bdytype;
//		int x_finest =  
		if (boundaries.south.bound(x)  && level == solver.L) {
			boundary_value = boundaries.south.inlet;
		}
		else {
			boundary_type = CLOSED;
		}
	}
	else if (nghbr_act_idx == -4 ) {

		sign = -1;
		dir = 2;
		edge = 4;
		boundary_type = boundaries.west.bdytype;
		if (boundaries.west.bound(y) && level == solver.L) {
			boundary_value = boundaries.west.inlet; // no need for q_source as we need discharge not h
		}
		else {
			boundary_type = CLOSED;
		}
	}



	h_owner = d_non_uniform_nghbrs.h_owner[idx]; // h0
	h_nghbr = d_non_uniform_nghbrs.h_nghbr[idx]; // h1
	z_owner = d_non_uniform_nghbrs.z_owner[idx]; // z0
	z_nghbr = d_non_uniform_nghbrs.z_nghbr[idx]; // z1


	//if (x == 40 && y == 1  && level == 7 /* && nghbr_act_idx == -1*/)
	//{
	//	printf
	//	(
	//		"h_nghbr: %f\n"
	//		"h_owner: %f\n"
	//		"nghbr_act_idx: %d\n"
	//		"nghbr_elem_idcs: %d\n"
	//		"boundary_type: %d\n"
	//		"boundary_value: %f\n"
	//		"edge: %d\n"
	//		"idx: %d\n\n"

	//		
	//		,
	//		h_nghbr, h_owner, nghbr_act_idx, d_non_uniform_nghbrs.nghbr_elem_idcs[idx], boundary_type, boundary_value, edge, idx
	//	);
	//}


	if (non_uniform_n) {
		NUMERIC_TYPE n_nghbr, n_owner;
		n_owner = d_non_uniform_nghbrs.n_owner[idx]; // n0
		n_nghbr = d_non_uniform_nghbrs.n_nghbr[idx]; // n1

		if (nghbr_act_idx < 0) {
			Manning = n_owner;
		}
		else {
			Manning = C(0.5) * (n_owner + n_nghbr);
		}

	}
	else {
		Manning = pars.FPn;
	}



	q_old = d_non_uniform_nghbrs.q[idx] / d_non_uniform_nghbrs.dx[idx];

	if (nghbr_act_idx < 0)
	{
		d_non_uniform_nghbrs.q[idx] = C(0.0); // default -> CLOSED
		
		if (boundary_type == FREE) 
		{
			hflow = h_owner;
			Sf = boundary_value; // 0.02; //BCptr->BC_Val[BCi]; // MKS must check this

			if (hflow > solver.DepthThresh)
			{
				d_non_uniform_nghbrs.q[idx] = sign * (FABS(q_old) + FABS(g * solver.Tstep * hflow * Sf)) / (C(1.0) + g * solver.Tstep * hflow * Manning * Manning * FABS(q_old) / (POW(hflow, (C(10.0) / C(3.0))))) * length;

				if (dir == 1 && sign == -1 && d_non_uniform_nghbrs.q[idx] > C(0.0)) d_non_uniform_nghbrs.q[idx] = C(0.0);
				if (dir == 1 && sign == 1 && d_non_uniform_nghbrs.q[idx] < C(0.0)) d_non_uniform_nghbrs.q[idx] = C(0.0);
				if (dir == 2 && sign == -1 && d_non_uniform_nghbrs.q[idx] > C(0.0)) d_non_uniform_nghbrs.q[idx] = C(0.0);
				if (dir == 2 && sign == 1 && d_non_uniform_nghbrs.q[idx] < C(0.0)) d_non_uniform_nghbrs.q[idx] = C(0.0);


			}
		}
		else if (boundary_type == HFIX || boundary_type == HVAR)
		{
			hflow = getmax(h_owner, boundary_value - z_owner);
			h0 = boundary_value;
			h1 = h_owner;
			z1 = z_owner;
			
			if (hflow >= solver.DepthThresh)
			{
				dh = h0 - z1 - h1;
				// change slops direction depending on the edge
				if (edge == 1 || edge == 4) {
					Sf = -dh / length;
				}
				else {
					Sf = dh / length;
				}
				// implement momentum equation
				d_non_uniform_nghbrs.q[idx] = (q_old - g * solver.Tstep * hflow * Sf) / (C(1.0) + g * solver.Tstep * hflow * Manning * Manning * FABS(q_old) / (POW(hflow, (C(10.0) / C(3.0))))) * length;
			}
		}
		else if (boundary_type == QFIX || boundary_type == QVAR)
		{
			d_non_uniform_nghbrs.q[idx] = -boundary_value * sign * length; // no need for q_source as we need discharge not h
		}

		if (d_non_uniform_nghbrs.q[idx] != 0)
		{
			d_non_uniform_nghbrs.v[idx] = d_non_uniform_nghbrs.q[idx] / length / hflow;
		}
		else d_non_uniform_nghbrs.v[idx] = C(0.0);


	}
	else 
	{
		if (h_owner > solver.DepthThresh || h_nghbr > solver.DepthThresh) {

			q_opp = C(0.0);

			int index = d_non_uniform_nghbrs.cumu_itrface_counts[idx];

			for (int i = 0; i < d_non_uniform_nghbrs.itrface_counts[idx]; i++)
			{
				q_opp += d_non_uniform_itrfaces.q_vol[index + i];
			}

			q_opp_avg = q_opp / d_non_uniform_nghbrs.dx_sum[idx];

			qvect = SQRT(q_old * q_old + q_opp_avg * q_opp_avg);

			y_owner = (z_owner + h_owner) * d_non_uniform_nghbrs.owner_flags[idx];
			y_nghbr = (z_nghbr + h_nghbr) * d_non_uniform_nghbrs.nghbr_flags[idx];

			dh = y_owner + y_nghbr; // y0 - y1

//			dx_owner = pars.dx * (1 << (solver.L - d_assem_sol.levels[d_non_uniform_nghbrs.owner_elem_idcs[idx]]));
//			dx_nghbr = pars.dx * (1 << (solver.L - d_assem_sol.levels[d_non_uniform_nghbrs.nghbr_elem_idcs[idx]]));

			Sf = -dh / length; // (C(0.5)* (dx_owner + dx_nghbr));
//			Sf = -dh / (C(0.5)* (dx_owner + dx_nghbr));
			hflow = FMAX(z_owner + h_owner, z_nghbr + h_nghbr) - FMAX(z_owner, z_nghbr);
			hflow = FMAX(hflow, C(0.0));
			MaxHflow = C(10.0);
			hflow = FMIN(hflow, MaxHflow);

			if (hflow > solver.DepthThresh)
			{
				d_non_uniform_nghbrs.q[idx] = ((theta * q_old) - (g * solver.Tstep * hflow * Sf)) / (C(1.0) + g * solver.Tstep * hflow * Manning * Manning * FABS(qvect) / (POW(hflow, (C(10.0) / C(3.0))))) * length; // Volumetric discharge Q
//				d_non_uniform_nghbrs.v = d_non_uniform_nghbrs.q[idx] / length / hflow;
			}
			else
			{
				d_non_uniform_nghbrs.q[idx] = C(0.0);
//				d_non_uniform_nghbrs.v = C(0.0);
			}
		}
		else
		{
			d_non_uniform_nghbrs.q[idx] = C(0.0);
//			d_non_uniform_nghbrs.v = C(0.0);
		}


		if (d_non_uniform_nghbrs.q[idx] != C(0.0))
		{
			d_non_uniform_nghbrs.v[idx] = d_non_uniform_nghbrs.q[idx] / length / hflow;
		}
		else d_non_uniform_nghbrs.v[idx] = C(0.0);



	}
}