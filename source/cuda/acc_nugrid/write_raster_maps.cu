#include "write_raster_maps.cuh"

__host__ void lis::cuda::acc_nugrid::write_raster_maps
(
	const char*                 respath,
	const AssembledSolution&    d_assem_sol,
	const int&                  mesh_dim,
	const Pars& pars,
	const int& call_gzip,
	const int& precision
)
{
	char fullpath[800];
	char tmp_sys_com[255];

	int nrows = 0;
	int ncols = 0;

	NUMERIC_TYPE cell_size = pars.dx;

	NUMERIC_TYPE xllcorner = pars.blx;
	NUMERIC_TYPE yllcorner = pars.bly;

	int NODATA_value = -9999;

	FILE* fp;
		
	ncols = pars.xsz;
	nrows = pars.ysz;

	NUMERIC_TYPE* h  = new NUMERIC_TYPE[mesh_dim * mesh_dim];
	NUMERIC_TYPE* z0  = new NUMERIC_TYPE[mesh_dim * mesh_dim];

	size_t bytes = mesh_dim * mesh_dim * sizeof(NUMERIC_TYPE);
	
	copy_cuda(h,	 d_assem_sol.h, bytes);
	copy_cuda(z0,  d_assem_sol.z0, bytes);

	// WRITING WATER DEPTH //

	snprintf(fullpath, 800 * sizeof(char), "%s-%.4d%s", respath, pars.SaveNo, ".wd");
	//
	fp = fopen(fullpath, "wb");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening file: lisflood.wd\n");
		exit(-1);
	}

	fprintf(fp, "ncols        %d\n", ncols);
	fprintf(fp, "nrows        %d\n", nrows);
	fprintf(fp, "xllcorner    %" NUM_FMT "\n", xllcorner);
	fprintf(fp, "yllcorner    %" NUM_FMT "\n", yllcorner);
	fprintf(fp, "cellsize     %" NUM_FMT "\n", cell_size);
	fprintf(fp, "NODATA_value %d\n", NODATA_value);
	
	for (int j = 0; j < nrows; j++)
	{
		for (int i = 0; i < ncols; i++)
		{
			int idx = (nrows - 1 - j) * mesh_dim + i;
//			int idx = j * mesh_dim + i;

			fprintf(fp, "%.*" NUM_FMT "\t", precision, (h[idx] < C(0.0)) ? C(0.0) : h[idx]);
		}

		fprintf(fp, "\n");
	}

	fclose(fp);

	// ------------------- //

	// WRITING WATER ELEV //

	//sprintf(fullpath, "%s%s%d%s", respath, "results-", saveint.count - 1, ".elev");

	//fp = fopen(fullpath, "w");

	//if (NULL == fp)
	//{
	//	fprintf(stderr, "Error opening file: lisflood.elev\n");
	//	exit(-1);
	//}

	//fprintf(fp, "ncols        %d\n", ncols);
	//fprintf(fp, "nrows        %d\n", nrows);
	//fprintf(fp, "xllcorner    %" NUM_FRMT "\n", xllcorner);
	//fprintf(fp, "yllcorner    %" NUM_FRMT "\n", yllcorner);
	//fprintf(fp, "cellsize     %" NUM_FRMT "\n", cell_size);
	//fprintf(fp, "NODATA_value %d\n", NODATA_value);

	//for (int j = 0; j < nrows; j++)
	//{
	//	for (int i = 0; i < ncols; i++)
	//	{
	//		int idx = j * mesh_dim + i;

	//		fprintf(fp, "%.15" NUM_FRMT " ", h[idx] + z[idx]);
	//	}

	//	fprintf(fp, "\n");
	//}

	//fclose(fp);

	// ------------------ //

	// WRITING X VELOCITIES //

	// reset the string https://stackoverflow.com/questions/1559487/how-to-empty-a-char-array
	//fullpath[0] = '\0';

	//sprintf(fullpath, "%s%s%d%s", respath, "results-", saveint.count - 1, ".Vx");

	//fp = fopen(fullpath, "w");

	//if (NULL == fp)
	//{
	//	fprintf(stderr, "Error opening file: lisflood.Vx\n");
	//	exit(-1);
	//}

	//fprintf(fp, "ncols        %d\n", ncols);
	//fprintf(fp, "nrows        %d\n", nrows);
	//fprintf(fp, "xllcorner    %" NUM_FRMT "\n", xllcorner);
	//fprintf(fp, "yllcorner    %" NUM_FRMT "\n", yllcorner);
	//fprintf(fp, "cellsize     %" NUM_FRMT "\n", cell_size);
	//fprintf(fp, "NODATA_value %d\n", NODATA_value);

	//for (int j = 0; j < nrows; j++)
	//{
	//	for (int i = 0; i < ncols; i++)
	//	{
	//		int idx = j * mesh_dim + i;

	//		NUMERIC_TYPE vx = qx[idx] / h[idx];

	//		fprintf(fp, "%.15" NUM_FRMT " ", vx);
	//	}

	//	fprintf(fp, "\n");
	//}

	//fclose(fp);

	// -------------------- //

	// WRITING Y VELOCITIES //

	// reset the string https://stackoverflow.com/questions/1559487/how-to-empty-a-char-array
	//fullpath[0] = '\0';

	//sprintf(fullpath, "%s%s%d%s", respath, "results-", saveint.count - 1, ".Vy");

	//fp = fopen(fullpath, "w");

	//if (NULL == fp)
	//{
	//	fprintf(stderr, "Error opening file: lisflood.Vy\n");
	//	exit(-1);
	//}

	//fprintf(fp, "ncols        %d\n", ncols);
	//fprintf(fp, "nrows        %d\n", nrows);
	//fprintf(fp, "xllcorner    %" NUM_FRMT "\n", xllcorner);
	//fprintf(fp, "yllcorner    %" NUM_FRMT "\n", yllcorner);
	//fprintf(fp, "cellsize     %" NUM_FRMT "\n", cell_size);
	//fprintf(fp, "NODATA_value %d\n", NODATA_value);

	//for (int j = 0; j < nrows; j++)
	//{
	//	for (int i = 0; i < ncols; i++)
	//	{
	//		int idx = j * mesh_dim + i;

	//		NUMERIC_TYPE vy = qy[idx] / h[idx];

	//		fprintf(fp, "%.15" NUM_FRMT " ", vy);
	//	}

	//	fprintf(fp, "\n");
	//}

	//fclose(fp);

	// -------------------- //

	// check if we need to zip the file up
	if (call_gzip == ON)
	{
		sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", fullpath);
		system(tmp_sys_com);
	}

	delete[] h;
	//delete[] qx;
	//delete[] qy;

	// ----------- //
}