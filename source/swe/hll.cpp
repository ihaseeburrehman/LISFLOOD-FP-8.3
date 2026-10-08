#include "hll.h"
#include "flux.h"
#include <cmath>
#include <algorithm>

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
)
{
	HLL(Solverptr, H_neg, HU_neg, HV_neg, H_pos, HU_pos, HV_pos,
			H_flux, HU_flux, HV_flux);
}

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
)
{
	HLL(Solverptr, H_neg, HV_neg, -HU_neg, H_pos, HV_pos, -HU_pos,
			H_flux, HV_flux, HU_flux);
	HU_flux = -HU_flux;
}

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
)
{
	NUMERIC_TYPE DepthThresh = Solverptr->DepthThresh;
	NUMERIC_TYPE g = Solverptr->g;

#include "hll_include.h"
}
