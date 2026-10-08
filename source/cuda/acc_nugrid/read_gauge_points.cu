#include "read_gauge_points.cuh"

void lis::cuda::acc_nugrid::read_gauge_points
(
	const char* stagefilename,
	const Pars& pars,
	const GaugePoints& gauge_points
)
{
	int num_gauge_points = 0;
	FILE* fp = fopen(stagefilename, "r");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening stage file for reading number of stage points, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}

	fscanf(fp, "%d", &num_gauge_points);

	NUMERIC_TYPE x_stage = C(0.0);
	NUMERIC_TYPE y_stage = C(0.0);

	for (int i = 0; i < gauge_points.num_points; i++)
	{
		fscanf(fp, "%" NUM_FMT " %" NUM_FMT, &x_stage, &y_stage);

//		int x =                  (x_stage - params.xmin) / params.cell_size;
//		int y = params.ysz - (y_stage - params.ymin) / params.cell_size;

		int x = (x_stage - pars.blx) / pars.dx;
		int y = (y_stage - pars.bly) / pars.dx;
		gauge_points.codes[i] = generate_morton_code(x, y);
	}

	fclose(fp);
}