#include "write_reals_to_file.cuh"

__host__ void lis::cuda::acc_nugrid::write_reals_to_file
(
	const char* filename,
	const char* respath,
	NUMERIC_TYPE*       d_results,
	const int&  array_length
)
{
	// allocating host array to copy_cuda device array to 
	NUMERIC_TYPE* h_results = new NUMERIC_TYPE[array_length];

	size_t bytes = array_length * sizeof(NUMERIC_TYPE);

	copy_cuda
	(
		h_results,
		d_results,
		bytes
	);

	FILE* fp;

	char fullpath[255];

	sprintf(fullpath, "%s%s%s", respath, filename, ".csv");

	fp = fopen(fullpath, "w");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening file : %s", filename);
		exit(-1);
	}

	fprintf(fp, "results\n");

	NUMERIC_TYPE sum = 0;

	for (int i = 0; i < array_length; i++)
	{
		fprintf(fp, "%.15" NUM_FMT "\n", h_results[i]);

		sum += h_results[i];
	}

	printf("%s: %.15f\n", filename, sum / array_length);

	fclose(fp);

	delete[] h_results;
}