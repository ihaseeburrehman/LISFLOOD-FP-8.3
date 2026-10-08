#include "write_mass_data.cuh"

void lis::cuda::acc_nugrid::write_mass_data
(
	const char* respath,
//	AssembledSolution& d_assem_sol,
	const GaugePoints& gauge_points,
	const int& mesh_dim,
	const Solver& solver,
	const States& states,
	const NUMERIC_TYPE& FArea,
	const NUMERIC_TYPE& vol2,
	const NUMERIC_TYPE& Qin,
	const NUMERIC_TYPE& Hds,
	const NUMERIC_TYPE& Qout,
	const NUMERIC_TYPE& Qerror,
	const NUMERIC_TYPE& Verror,
	const NUMERIC_TYPE& RainTotalLoss,
	const NUMERIC_TYPE& InfilTotalLoss,
	const NUMERIC_TYPE& EvapTotalLoss
)
{
//	if (gauge_points.num_points == 0) return;

//	int num_finest_elems = mesh_dim * mesh_dim;

//	size_t bytes = sizeof(NUMERIC_TYPE) * num_finest_elems;

//	NUMERIC_TYPE* h = new NUMERIC_TYPE[num_finest_elems];
//	NUMERIC_TYPE* qx = new NUMERIC_TYPE[num_finest_elems];
//	NUMERIC_TYPE* qy = new NUMERIC_TYPE[num_finest_elems];

//	copy_cuda(h, d_assem_sol.h, sizeof(NUMERIC_TYPE) * num_finest_elems);
//	copy_cuda(qx, d_assem_sol.qx0, sizeof(NUMERIC_TYPE) * num_finest_elems);
//	copy_cuda(qy, d_assem_sol.qy0, sizeof(NUMERIC_TYPE) * num_finest_elems);

	char fullpath[255];

	sprintf(fullpath, "%s%s", respath, ".mass");

	FILE* fp; // = (solver.t > 0) ? fopen(fullpath, "a") : fopen(fullpath, "w");

	if (solver.t > 0) { //if this is a checkpointed job, we only need to amend the .mass file
		fp = fopen(fullpath, "a");
	}
	else {
		fp = fopen(fullpath, "w");
	}

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening mass balance file, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}

	if (solver.t == C(0.0))
	{
		fprintf(fp, "Time         Tstep      MinTstep   NumTsteps    Area         Vol         Qin         Hds        Qout          Qerror       Verror       Rain-(Inf+Evap)\n");
	}
	else 
	{

		fprintf(fp, "%-12.3" NUM_FMT" %-10.4" NUM_FMT" %-10.4" NUM_FMT" %-10ld %12.4" NUM_FMT" %12.4" NUM_FMT"  %-11.3" NUM_FMT" %-10.3" NUM_FMT" %-11.3" NUM_FMT" %12.4" NUM_FMT" %12.4" NUM_FMT" %12.4" NUM_FMT"\n", solver.t, solver.Tstep, solver.MinTstep, solver.itCount, FArea, vol2, Qin, Hds, Qout, Qerror, Verror, RainTotalLoss - (InfilTotalLoss + EvapTotalLoss));

	}

	fclose(fp);

//	delete[] h;
//	delete[] qx;
//	delete[] qy;
}