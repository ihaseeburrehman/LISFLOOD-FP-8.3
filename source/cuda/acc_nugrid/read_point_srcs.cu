#include "read_point_srcs.cuh"

void lis::cuda::acc_nugrid::read_point_srcs
(
	const char* bcifilename,
	const char* bdyfilename,
	const Pars& pars,
	const NUMERIC_TYPE& time_now,
	PointSources& point_sources,
	const BoundCs& boundCs
)
{
	char str[255];

	if (strlen(bcifilename) == 0)
	{
		//		printf("No point sources counted in boundary condition file, proceeding with zero point sources.\n");
		return;
	}

	FILE* fp = fopen(bcifilename, "r");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening point source file, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}



	NUMERIC_TYPE x_stage = C(0.0);
	NUMERIC_TYPE y_stage = C(0.0);

	int x = 0;
	int y = 0;

	char point = '\0';
	char inlet_type_buf[8] = { '\0' };
	int  srcs_counted = 0;

	while (srcs_counted < point_sources.num_srcs)
	{
		fgets(str, sizeof(str), fp);

		sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s", &point, &x_stage, &y_stage, inlet_type_buf);

		if (point == 'P')
		{
//			x =                  (x_stage - params.xmin) / params.cell_size;
//			y = params.ysz - (y_stage - params.ymin) / params.cell_size;
			x = (x_stage - pars.blx) / pars.dx;
			y = (y_stage - pars.bly) / pars.dx;

			point_sources.h_codes[srcs_counted] = generate_morton_code(x, y);

			if (!strncmp(inlet_type_buf, "HFIX", 4))
			{
				point_sources.h_src_types[srcs_counted] = HFIX;
				sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %" NUM_FMT, &point, &x_stage, &y_stage, inlet_type_buf, &point_sources.h_srcs[srcs_counted]);
			}
			else if (!strncmp(inlet_type_buf, "HVAR", 4))
			{
				point_sources.h_src_types[srcs_counted] = HVAR;
				sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %s", &point, &x_stage, &y_stage, inlet_type_buf, &point_sources.timeseries[srcs_counted * 32]);
			}
			else if (!strncmp(inlet_type_buf, "QFIX", 4))
			{
				point_sources.h_src_types[srcs_counted] = QFIX;
				sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %" NUM_FMT, &point, &x_stage, &y_stage, inlet_type_buf, &point_sources.h_srcs[srcs_counted]);
			}
			else if (!strncmp(inlet_type_buf, "QVAR", 4))
			{
				point_sources.h_src_types[srcs_counted] = QVAR;
				sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %s", &point, &x_stage, &y_stage, inlet_type_buf, &point_sources.timeseries[srcs_counted * 32]);
			}

			srcs_counted++;
		}
	}

	fclose(fp);

	size_t bytes_codes = sizeof(MortonCode) * point_sources.num_srcs;
	size_t bytes_srcs = sizeof(NUMERIC_TYPE) * point_sources.num_srcs;
	size_t bytes_src_types = sizeof(int) * point_sources.num_srcs;

	copy_cuda(point_sources.d_codes, point_sources.h_codes, bytes_codes);
	copy_cuda(point_sources.d_srcs, point_sources.h_srcs, bytes_srcs);
	copy_cuda(point_sources.d_src_types, point_sources.h_src_types, bytes_src_types);

	point_sources.update_all_sources(bdyfilename, time_now, boundCs);
}