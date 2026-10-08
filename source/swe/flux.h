#pragma once
#include "../lisflood.h"

void physical_flux_x
(
	Solver *Solverptr,
	NUMERIC_TYPE H,
	NUMERIC_TYPE HU,
	NUMERIC_TYPE HV,
	NUMERIC_TYPE& H_flux,
	NUMERIC_TYPE& HU_flux,
	NUMERIC_TYPE& HV_flux
);

void physical_flux_y
(
	Solver *Solverptr,
	NUMERIC_TYPE H,
	NUMERIC_TYPE HU,
	NUMERIC_TYPE HV,
	NUMERIC_TYPE& H_flux,
	NUMERIC_TYPE& HU_flux,
	NUMERIC_TYPE& HV_flux
);

void physical_flux
(
	Solver *Solverptr,
	NUMERIC_TYPE H,
	NUMERIC_TYPE HU,
	NUMERIC_TYPE HV,
	NUMERIC_TYPE& H_flux,
	NUMERIC_TYPE& HU_flux,
	NUMERIC_TYPE& HV_flux
);
