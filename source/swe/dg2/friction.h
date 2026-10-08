#pragma once
#include "../../lisflood.h"

namespace dg2
{
	void apply_friction
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr
	);

	void friction
	(
		Solver *Solverptr,
		NUMERIC_TYPE H,
		NUMERIC_TYPE& HU,
		NUMERIC_TYPE& HV,
        NUMERIC_TYPE n
	);
}
