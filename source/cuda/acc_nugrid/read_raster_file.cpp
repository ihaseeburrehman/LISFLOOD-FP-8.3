#include "read_raster_file.h"
#include <cmath>

void lis::cuda::acc_nugrid::read_raster_file
(
	const char* raster_filename,
	NUMERIC_TYPE* raster_array,
	const int&  mesh_dim,
	const NUMERIC_TYPE& no_data
)
{
	char buf[255];

	int nrows = 0;
	int ncols = 0;

	NUMERIC_TYPE dummy = C(0.0);
	NUMERIC_TYPE nodata = C(0.0);

	FILE* fp = fopen(raster_filename, "r");

	if (NULL == fp)
	{
		fprintf(stdout, "No raster file found: %s, using default values.\n", raster_filename);
		return;
	}

	fscanf(fp, "%s %d", buf, &ncols);
	fscanf(fp, "%s %d", buf, &nrows);
	fscanf(fp, "%s %" NUM_FMT, buf, &dummy);
	fscanf(fp, "%s %" NUM_FMT, buf, &dummy);
	fscanf(fp, "%s %" NUM_FMT, buf, &dummy);
	fscanf(fp, "%s %" NUM_FMT, buf, &nodata);

	for (int j = 0; j < nrows; j++)
	{
		for (int i = 0; i < ncols; i++)
		{
			fscanf( fp, "%" NUM_FMT, &raster_array[(nrows - 1 - j) * mesh_dim + i] );
//			fscanf( fp, "%" NUM_FMT, &raster_array[j * mesh_dim + i] );
			if (FABS(raster_array[(nrows - 1 - j) * mesh_dim + i] - nodata) < 1e-5 )
			{
				raster_array[(nrows - 1 - j) * mesh_dim + i] = no_data;
			}

		}
	}

	fclose(fp);
}