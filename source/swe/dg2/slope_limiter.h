#pragma once
#include "../../lisflood.h"
#include "../dg2.h"

namespace dg2
{
    void zero_perimeter_slopes
    (
		Pars *Parptr,
        Arrays *Arrptr,
		FlowCoefficients const& U
    );

	void apply_slope_limiter
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U
	);

	void apply_slope_limiter_x
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U
	);

	void apply_slope_limiter_y
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U
	);

	NUMERIC_TYPE limit_slope
	(
		Solver *Solverptr,
		NUMERIC_TYPE stencil_minH,
		NUMERIC_TYPE mesh_delta,
		NUMERIC_TYPE backward_const,
		NUMERIC_TYPE backward_slope,
		NUMERIC_TYPE local_const,
		NUMERIC_TYPE local_slope,
		NUMERIC_TYPE forward_const,
		NUMERIC_TYPE forward_slope
	);

	NUMERIC_TYPE minmod(NUMERIC_TYPE a, NUMERIC_TYPE b, NUMERIC_TYPE c);

	NUMERIC_TYPE stencil_minH_x
	(
		Pars *Parptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE stencil_minH_y
	(
		Pars *Parptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);
}
