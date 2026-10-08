#pragma once

#include "cuda_runtime.h" //TOREMOVE
#include "Boundary.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct Boundaries
{
	Boundary north;
	Boundary east;
	Boundary south;
	Boundary west;

	Boundaries
	(
		const Fnames& filenames,
		const Pars& pars
	)
	:
		north(filenames, pars, NORTH),
		east (filenames, pars, EAST),
		south(filenames, pars, SOUTH),
		west (filenames, pars, WEST)
	{
//		if (test_case != 0) fprintf(stdout, "Running built-in test case, using open boundary conditions.\n");
	}

	void update_all_inlets
	(
		const char* input_filename,
		const NUMERIC_TYPE& time_now
	)
	{
		north.update_inlet(time_now);
		east.update_inlet (time_now);
		south.update_inlet(time_now);
		west.update_inlet (time_now);
	}

} Boundaries;

}
}
}