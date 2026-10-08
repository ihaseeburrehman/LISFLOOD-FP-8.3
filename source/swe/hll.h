#pragma once
#include "../lisflood.h"

void HLL_x
(
	Solver *Solverptr,
	NUMERIC_TYPE H_neg,
	NUMERIC_TYPE HU_neg,
	NUMERIC_TYPE HV_neg,
	NUMERIC_TYPE H_pos,
	NUMERIC_TYPE HU_pos,
	NUMERIC_TYPE HV_pos,
	NUMERIC_TYPE& H_flux,
	NUMERIC_TYPE& HU_flux,
	NUMERIC_TYPE& HV_flux
);

void HLL_y
(
	Solver *Solverptr,
	NUMERIC_TYPE H_neg,
	NUMERIC_TYPE HU_neg,
	NUMERIC_TYPE HV_neg,
	NUMERIC_TYPE H_pos,
	NUMERIC_TYPE HU_pos,
	NUMERIC_TYPE HV_pos,
	NUMERIC_TYPE& H_flux,
	NUMERIC_TYPE& HU_flux,
	NUMERIC_TYPE& HV_flux
);

void HLL
(
	Solver *Solverptr,
	NUMERIC_TYPE H_neg,
	NUMERIC_TYPE HU_neg,
	NUMERIC_TYPE HV_neg,
	NUMERIC_TYPE H_pos,
	NUMERIC_TYPE HU_pos,
	NUMERIC_TYPE HV_pos,
	NUMERIC_TYPE& H_flux,
	NUMERIC_TYPE& HU_flux,
	NUMERIC_TYPE& HV_flux
);
