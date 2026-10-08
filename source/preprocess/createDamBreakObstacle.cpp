#include "../lisflood.h"
#include "../utility.h"
#include "../swe/dg2/fields.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <initializer_list>

NUMERIC_TYPE Z(NUMERIC_TYPE x, NUMERIC_TYPE y)
{
	NUMERIC_TYPE peak = C(0.4);
	if (x > C(25.5) && x <= C(28.5))
	{
		return peak/C(3.0) * (x - C(25.5));
	}
	else if (x > C(28.5) && x < C(31.5))
	{
		return -peak/C(3.0) * (x - C(31.5));
	}
	else
	{
		return C(0.0);
	}
}

NUMERIC_TYPE H(NUMERIC_TYPE x, NUMERIC_TYPE y)
{
	if (x < C(15.5))
	{
		return C(0.75);
	}
	else
	{
		return C(0.0);
	}
}

int main()
{
	Pars pars;
	pars.xsz = 114;
	pars.ysz = 4;
	pars.blx = C(0.0);
	pars.bly = C(0.0);
	pars.dx = C(1.0)/C(3.0);
	pars.dy = pars.dx;
	Pars *Parptr = &pars;

	Arrays arrays;
	Arrays *Arrptr = &arrays;

	Arrptr->DEM = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);
	Arrptr->H = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);
	dg2::allocate_fields(Parptr, Arrptr);

	dg2::initialise_field(Z, Parptr, Arrptr->DEM, Arrptr->DEM1x, Arrptr->DEM1y);
	dg2::initialise_field(H, Parptr, Arrptr->H, Arrptr->H1x, Arrptr->H1y);

	dg2::write_dem("damBreakObstacle.dem", Parptr, Arrptr);
	dg2::write_startfile("damBreakObstacle.start", Parptr, Arrptr);

	dg2::deallocate_fields(Parptr, Arrptr);
	return EXIT_SUCCESS;
}
