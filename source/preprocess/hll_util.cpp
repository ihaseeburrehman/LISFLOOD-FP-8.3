#include <cstdlib>
#include "../swe/hll.h"

int main(int argc, char *argv[])
{
	if (argc != 7)
	{
		fprintf(stderr, "Usage: %s H_neg HU_neg HV_neg H_pos HU_pos HV_pos\n",
				argv[0]);
		return EXIT_FAILURE;
	}

	Solver solver;
	solver.g = C(9.80665);
	solver.DepthThresh = C(1e-3);

	NUMERIC_TYPE H_neg = atof(argv[1]);
	NUMERIC_TYPE HU_neg = atof(argv[2]);
	NUMERIC_TYPE HV_neg = atof(argv[3]);
	NUMERIC_TYPE H_pos = atof(argv[4]);
	NUMERIC_TYPE HU_pos = atof(argv[5]);
	NUMERIC_TYPE HV_pos = atof(argv[6]);

	NUMERIC_TYPE FH = C(0.0);
	NUMERIC_TYPE FHU = C(0.0);
	NUMERIC_TYPE FHV = C(0.0);

	HLL_x(&solver, H_neg, HU_neg, HV_neg, H_pos, HU_pos, HV_pos, FH, FHU, FHV);

	printf("%.17g %.17g %.17g\n", FH, FHU, FHV);

	return EXIT_SUCCESS;
}

