#include "read_num_stage_points.h"

int lis::cuda::acc_nugrid::read_num_stage_points(const char* stagefilename)
{
	int num_stage_points = 0;

	FILE*  fp = fopen(stagefilename, "r");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening stage file for reading number of stage points, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}

	fscanf(fp, "%d", &num_stage_points);

	fclose(fp);

	return num_stage_points;
}