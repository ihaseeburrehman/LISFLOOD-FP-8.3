#include "../lisflood.h"
#include "../utility.h"
#include "../swe/dg2/fields.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <initializer_list>

NUMERIC_TYPE Z(NUMERIC_TYPE x, NUMERIC_TYPE y)
{
	if (C(11.0) <= y && y <= C(19.0))
	{
		if (C(16.0) <= x && x <= C(24.0))
		{
			return C(0.86);
		}
		else if (C(36.0) <= x && x <= C(44.0))
		{
			return C(1.78);
		}
		else if (C(56.0) <= x && x <= C(64.0))
		{
			return C(2.3);
		}
	}

	return C(0.0);
}

int main()
{
	Pars pars;
	pars.xsz = 75;
	pars.ysz = 30;
	pars.blx = C(0.0);
	pars.bly = C(0.0);
	pars.tly = C(30.0);
	pars.dx = C(1.0);
	pars.dy = pars.dx;
	Pars *Parptr = &pars;

	Arrays arrays;
	Arrays *Arrptr = &arrays;

	Arrptr->DEM = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);
	Arrptr->H = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);
	dg2::allocate_fields(Parptr, Arrptr);

	dg2::initialise_field(Z, Parptr, Arrptr->DEM, Arrptr->DEM1x, Arrptr->DEM1y);
	dg2::initialise_h_from_eta(C(1.95), Parptr, Arrptr);

	dg2::write_dem("lakeAtRestBlocks.dem", Parptr, Arrptr);
	dg2::write_startfile("lakeAtRestBlocks.start", Parptr, Arrptr);

	dg2::deallocate_fields(Parptr, Arrptr);
	return EXIT_SUCCESS;
}
