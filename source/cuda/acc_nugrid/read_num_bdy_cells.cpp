#include "read_num_bdy_cells.h"

int lis::cuda::acc_nugrid::read_num_bdy_cells
(
	const char* bcifilename,
	const Pars& pars,
	const int                   direction
)
{
	NUMERIC_TYPE upper = C(0.0);
	NUMERIC_TYPE lower = C(0.0);

	char bcidir = '0';
	char filedir = '0';

	char str[255];
	char bdytype_buf[8];

	NUMERIC_TYPE origin = C(0.0);

	switch (direction)
	{
	case NORTH:
		origin = pars.blx;
		bcidir = 'N';
		break;
	case EAST:
		origin = pars.bly;
		bcidir = 'E';
		break;
	case SOUTH:
		origin = pars.blx;
		bcidir = 'S';
		break;
	case WEST:
		origin = pars.bly;
		bcidir = 'W';
		break;
	default:
		break;
	}

	FILE* fp = fopen(bcifilename, "r");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening boundary condition file, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}

	while (filedir != bcidir)
	{
		if (NULL == fgets(str, sizeof(str), fp))
		{
			fprintf(stdout, "No enforced boundary cells counted for boundary %c.\n", bcidir);
			fclose(fp);

			return 0;
		}

		sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s", &filedir, &lower, &upper, bdytype_buf);
	}

	fclose(fp);

	int start = (lower - origin) / pars.dx;
	int end = (upper - origin) / pars.dx - 1;

	return end - start + 1;
}