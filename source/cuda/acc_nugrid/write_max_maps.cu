#include "write_max_maps.cuh"

__host__ void lis::cuda::acc_nugrid::write_max_maps
(
	const char*                 respath,
	const int&                  mesh_dim,
	const Pars& pars,
	NUMERIC_TYPE* d_maxH,
	NUMERIC_TYPE* d_totalHtm,
	NUMERIC_TYPE* d_maxHtm,
	NUMERIC_TYPE* d_initHtm,
	const int&                  precision
)
{
	char fullpath[800];

	int nrows = 0;
	int ncols = 0;

	NUMERIC_TYPE cell_size = pars.dx;

	NUMERIC_TYPE xllcorner = pars.blx;
	NUMERIC_TYPE yllcorner = pars.bly;

	int NODATA_value = -9999;

	FILE* fp;
		
	ncols = pars.xsz;
	nrows = pars.ysz;

	NUMERIC_TYPE* maxH = new NUMERIC_TYPE[mesh_dim * mesh_dim];
	NUMERIC_TYPE* totalHtm = new NUMERIC_TYPE[mesh_dim * mesh_dim];
	NUMERIC_TYPE* maxHtm = new NUMERIC_TYPE[mesh_dim * mesh_dim];
	NUMERIC_TYPE* initHtm = new NUMERIC_TYPE[mesh_dim * mesh_dim];

	size_t bytes = mesh_dim * mesh_dim * sizeof(NUMERIC_TYPE);

	copy_cuda(maxH,     d_maxH,     bytes);
	copy_cuda(totalHtm, d_totalHtm, bytes);
	copy_cuda(maxHtm,   d_maxHtm,   bytes);
	copy_cuda(initHtm,  d_initHtm,  bytes);


	// WRITING .inittm //

	snprintf(fullpath, 800 * sizeof(char), "%s%s", respath, ".inittm");

	//
	fp = fopen(fullpath, "wb");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening file: .inittm\n");
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

			fprintf(fp, "%.*" NUM_FMT "\t", precision, initHtm[idx]);

		}

		fprintf(fp, "\n");
	}

	fclose(fp);

	// ------------------- //

	// WRITING .totaltm //

	snprintf(fullpath, 800 * sizeof(char), "%s%s", respath, ".totaltm");

	//
	fp = fopen(fullpath, "wb");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening file: .totaltm\n");
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

			fprintf(fp, "%.*" NUM_FMT "\t", precision, totalHtm[idx]);

		}

		fprintf(fp, "\n");
	}

	fclose(fp);

	// -------------------- //

	// WRITING .max //

	snprintf(fullpath, 800 * sizeof(char), "%s%s", respath, ".max");

	//
	fp = fopen(fullpath, "wb");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening file: .max\n");
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

			fprintf(fp, "%.*" NUM_FMT "\t", precision, maxH[idx]);

		}

		fprintf(fp, "\n");
	}

	fclose(fp);

	// -------------------- //

	// WRITING .maxtm //

	snprintf(fullpath, 800 * sizeof(char), "%s%s", respath, ".maxtm");

	//
	fp = fopen(fullpath, "wb");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening file: .maxtm\n");
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

			fprintf(fp, "%.*" NUM_FMT "\t", precision, maxHtm[idx]);

		}

		fprintf(fp, "\n");
	}

	fclose(fp);

	// -------------------- //


	delete[] maxH;
	delete[] totalHtm;
	delete[] maxHtm;
	delete[] initHtm;


	// ----------- //
}