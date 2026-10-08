#include "flux.h"

void physical_flux_x
(
	Solver *Solverptr,
	NUMERIC_TYPE H,
	NUMERIC_TYPE HU,
	NUMERIC_TYPE HV,
	NUMERIC_TYPE& H_flux,
	NUMERIC_TYPE& HU_flux,
	NUMERIC_TYPE& HV_flux
)
{
	physical_flux(Solverptr, H, HU, HV, H_flux, HU_flux, HV_flux);
}

void physical_flux_y
(
	Solver *Solverptr,
	NUMERIC_TYPE H,
	NUMERIC_TYPE HU,
	NUMERIC_TYPE HV,
	NUMERIC_TYPE& H_flux,
	NUMERIC_TYPE& HU_flux,
	NUMERIC_TYPE& HV_flux
)
{
	physical_flux(Solverptr, H, HV, HU, H_flux, HV_flux, HU_flux);
}

void physical_flux
(
	Solver *Solverptr,
	NUMERIC_TYPE H,
	NUMERIC_TYPE HU,
	NUMERIC_TYPE HV,
	NUMERIC_TYPE& H_flux,
	NUMERIC_TYPE& HU_flux,
	NUMERIC_TYPE& HV_flux
)
{
	if (H <= Solverptr->DepthThresh)
	{
		H_flux = C(0.0);
		HU_flux = C(0.0);
		HV_flux = C(0.0);
	}
	else
	{
		H_flux = HU;
		HU_flux = HU*HU/H + C(0.5)*Solverptr->g*H*H;
		HV_flux = HU*HV/H;
	}
}
