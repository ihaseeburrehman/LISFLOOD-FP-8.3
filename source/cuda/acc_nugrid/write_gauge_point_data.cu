#include "write_gauge_point_data.cuh"

void lis::cuda::acc_nugrid::write_gauge_point_data
(
	const char* respath,
	AssembledSolution& d_assem_sol,
	const GaugePoints& gauge_points,
	const int& mesh_dim,
	const char* stagefilename,
	const Pars& pars,
	const Solver& solver,
	const States& states
)
{
	if (gauge_points.num_points == 0) return;

	int num_finest_elems = mesh_dim * mesh_dim;
	
	size_t bytes = sizeof(NUMERIC_TYPE) * num_finest_elems;

	NUMERIC_TYPE* h = new NUMERIC_TYPE[num_finest_elems];
	NUMERIC_TYPE* z = new NUMERIC_TYPE[num_finest_elems];
//	NUMERIC_TYPE* qx = new NUMERIC_TYPE[num_finest_elems];
//	NUMERIC_TYPE* qy = new NUMERIC_TYPE[num_finest_elems];

	copy_cuda(h, d_assem_sol.h, bytes);
	copy_cuda(z, d_assem_sol.z0, bytes);
//	copy_cuda(qx, d_assem_sol.qx0, bytes);
//	copy_cuda(qy, d_assem_sol.qy0, bytes);
	
	char fullpath[255];

	sprintf(fullpath, "%s%s", respath, ".stage");

	FILE* fp; // = (time_now > 0) ? fopen(fullpath, "a") : fopen(fullpath, "w");

	if (/*states.checkpoint == ON && */ solver.t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
		fp = fopen(fullpath, "a");
	}
	else {
		fp = fopen(fullpath, "w");
	}

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening stage results file, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}

	if (solver.t == C(0.0) /* || states.checkpoint == OFF*/)
	{
		//fprintf(fp, "time,");

		//for (int i = 0; i < gauge_points.num_points; i++) fprintf(fp, "stage%d,", i + 1);

		//fprintf(fp, "\n");

		fprintf(fp, "Stage output, depth (m). Stage locations from: %s\n\n", stagefilename);
		fprintf(fp, "Stage information (stage,x,y,elev):\n");
		for (int i = 0; i < gauge_points.num_points; i++)
		{
			int x = compact(gauge_points.codes[i]);
			int y = compact(gauge_points.codes[i] >> 1);

			NUMERIC_TYPE x_centre = pars.blx + x * pars.dx + pars.dx / C(2.0);
			NUMERIC_TYPE y_centre = pars.bly + y * pars.dx + pars.dx / C(2.0);

			index_1D idx = y * mesh_dim + x;

			fprintf(fp, "%d\t%.4" NUM_FMT"\t%.4" NUM_FMT"\t%.4" NUM_FMT"\n", i + 1, x_centre, y_centre, z[idx]);
		}

		fprintf(fp, "\nOutput, depths:\n");
		fprintf(fp, "Time; stages 1 to %d\n", gauge_points.num_points);


	}
	else {

		//fprintf(fp, "%" NUM_FMT ",", time_now);

		//for (int i = 0; i < gauge_points.num_points; i++)
		//{
		//	int x = compact(gauge_points.codes[i]);
		//	int y = compact(gauge_points.codes[i] >> 1);

		//	index_1D idx = y * mesh_dim + x;

		//	fprintf(fp, "%" NUM_FMT ",", h[idx]);
		//}

		//fprintf(fp, "\n");

		fprintf(fp, "%12.3" NUM_FMT"", solver.t);
		for (int i = 0; i < gauge_points.num_points; i++)
		{

			int x = compact(gauge_points.codes[i]);
			int y = compact(gauge_points.codes[i] >> 1);

			index_1D idx = y * mesh_dim + x;

			fprintf(fp, "%10.4" NUM_FMT"", h[idx]);
		}
		fprintf(fp, "\n");

		//fflush(Fptr->stage_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
		//// added to export scalar velocity
		//if (Statesptr->voutput_stage == ON)
		//{
		//	fprintf(Fptr->vel_fp, "%12.3" NUM_FMT"", Solverptr->t);
		//	for (i = 0; i < Locptr->Nstages; i++)
		//	{
		//		int index = Locptr->stage_grid_x[i] + Locptr->stage_grid_y[i] * (Parptr->xsz + 1);
		//		if (Locptr->stage_check[i] == 1) fprintf(Fptr->vel_fp, "%10.4" NUM_FMT"", sqrt(pow(getmax(fabs(Arrptr->Vx[index]), fabs(Arrptr->Vx[index + 1])), 2) + pow(getmax(fabs(Arrptr->Vy[index]), fabs(Arrptr->Vy[index + (Parptr->xsz + 1)])), 2)));
		//		else fprintf(Fptr->vel_fp, "-\t");
		//	}
		//	fprintf(Fptr->vel_fp, "\n");
		//	fflush(Fptr->vel_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
		//}

	}
	fclose(fp);

	delete[] h;
	delete[] z;
//	delete[] qx;
//	delete[] qy;
}

void lis::cuda::acc_nugrid::write_velocity_point_data
(
	const char* respath,
	AssembledSolution& d_assem_sol,
	const GaugePoints& gauge_points,
	const int& mesh_dim,
	const char* stagefilename,
	const Pars& pars,
	const Solver& solver,
	const States& states
)
{
	if (gauge_points.num_points == 0) return;

	int num_finest_elems = mesh_dim * mesh_dim;

    size_t bytes = sizeof(NUMERIC_TYPE) * num_finest_elems;

	NUMERIC_TYPE* v = new NUMERIC_TYPE[num_finest_elems];
	NUMERIC_TYPE* z = new NUMERIC_TYPE[num_finest_elems];
	//	NUMERIC_TYPE* qx = new NUMERIC_TYPE[num_finest_elems];
	//	NUMERIC_TYPE* qy = new NUMERIC_TYPE[num_finest_elems];

	copy_cuda(v, d_assem_sol.v, bytes);
	copy_cuda(z, d_assem_sol.z0, bytes);
	//	copy_cuda(qx, d_assem_sol.qx0, bytes);
	//	copy_cuda(qy, d_assem_sol.qy0, bytes);

	char fullpath[255];

	sprintf(fullpath, "%s%s", respath, ".velocity");

	FILE* fp; // = (time_now > 0) ? fopen(fullpath, "a") : fopen(fullpath, "w");

	if (/*states.checkpoint == ON && */ solver.t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
		fp = fopen(fullpath, "a");
	}
	else {
		fp = fopen(fullpath, "w");
	}

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening stage results file, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}

	if (solver.t == C(0.0) /* || states.checkpoint == OFF*/)
	{
		//fprintf(fp, "time,");

		//for (int i = 0; i < gauge_points.num_points; i++) fprintf(fp, "stage%d,", i + 1);

		//fprintf(fp, "\n");

		fprintf(fp, "Velocity output, velocity (ms-1). Velocity locations from: %s\n\n", stagefilename);
		fprintf(fp, "Stage information (stage,x,y,elev):\n");
		for (int i = 0; i < gauge_points.num_points; i++)
		{
			int x = compact(gauge_points.codes[i]);
			int y = compact(gauge_points.codes[i] >> 1);

			NUMERIC_TYPE x_centre = pars.blx + x * pars.dx + pars.dx / C(2.0);
			NUMERIC_TYPE y_centre = pars.bly + y * pars.dx + pars.dx / C(2.0);

			index_1D idx = y * mesh_dim + x;

			fprintf(fp, "%d\t%.4" NUM_FMT"\t%.4" NUM_FMT"\t%.4" NUM_FMT"\n", i + 1, x_centre, y_centre, z[idx]);
		}

		fprintf(fp, "\nOutput, depths:\n");
		fprintf(fp, "Time; velocities 1 to %d\n", gauge_points.num_points);


	}
	else {

		//fprintf(fp, "%" NUM_FMT ",", time_now);

		//for (int i = 0; i < gauge_points.num_points; i++)
		//{
		//	int x = compact(gauge_points.codes[i]);
		//	int y = compact(gauge_points.codes[i] >> 1);

		//	index_1D idx = y * mesh_dim + x;

		//	fprintf(fp, "%" NUM_FMT ",", h[idx]);
		//}

		//fprintf(fp, "\n");

		fprintf(fp, "%12.3" NUM_FMT"", solver.t);
		for (int i = 0; i < gauge_points.num_points; i++)
		{

			int x = compact(gauge_points.codes[i]);
			int y = compact(gauge_points.codes[i] >> 1);

			index_1D idx = y * mesh_dim + x;

			fprintf(fp, "%10.4" NUM_FMT"", v[idx]);
		}
		fprintf(fp, "\n");

		//fflush(Fptr->stage_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
		//// added to export scalar velocity
		//if (Statesptr->voutput_stage == ON)
		//{
		//	fprintf(Fptr->vel_fp, "%12.3" NUM_FMT"", Solverptr->t);
		//	for (i = 0; i < Locptr->Nstages; i++)
		//	{
		//		int index = Locptr->stage_grid_x[i] + Locptr->stage_grid_y[i] * (Parptr->xsz + 1);
		//		if (Locptr->stage_check[i] == 1) fprintf(Fptr->vel_fp, "%10.4" NUM_FMT"", sqrt(pow(getmax(fabs(Arrptr->Vx[index]), fabs(Arrptr->Vx[index + 1])), 2) + pow(getmax(fabs(Arrptr->Vy[index]), fabs(Arrptr->Vy[index + (Parptr->xsz + 1)])), 2)));
		//		else fprintf(Fptr->vel_fp, "-\t");
		//	}
		//	fprintf(Fptr->vel_fp, "\n");
		//	fflush(Fptr->vel_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
		//}

	}
	fclose(fp);

	delete[] v;
	delete[] z;
	//	delete[] qx;
	//	delete[] qy;
}