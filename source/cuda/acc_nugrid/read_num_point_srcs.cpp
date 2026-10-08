#include "read_num_point_srcs.h"

lis::cuda::acc_nugrid::PointSources lis::cuda::acc_nugrid::read_num_point_srcs
(
	const char* bcifilename
)
{
	char str[255];

	if (strlen(bcifilename) == 0)
	{
		//		printf("No point sources counted in boundary condition file, proceeding with zero point sources.\n");
		return 0;
	}

	FILE* fp = fopen(bcifilename, "r");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening point source file for counting number of point sources, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}


	char point = '\0';
	int  num_srcs = 0;

	NUMERIC_TYPE x_stage = C(0.0); //TOREMOVE
	NUMERIC_TYPE y_stage = C(0.0); //TOREMOVE

	while (!(NULL == fgets(str, sizeof(str), fp)))
	{
		sscanf(str, "%c", &point);

		if (point == 'P') num_srcs++;
	}

	if (num_srcs == 0)
	{
		fprintf(stdout, "No point sources counted in boundary condition file, proceeding with zero point sources.\n");
		fclose(fp);

		return 0;
	}

	fclose(fp);

	return num_srcs;
}