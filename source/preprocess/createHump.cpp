#include "../lisflood.h"
#include "../utility.h"
#include "../swe/dg2/fields.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <initializer_list>

NUMERIC_TYPE Z(NUMERIC_TYPE x, NUMERIC_TYPE y)
{
	NUMERIC_TYPE a = C(0.6);
	NUMERIC_TYPE lambda = C(10.0);

	return a * pow(C(1.0) / cosh(M_PI*x/lambda), 2);
}

int main()
{
	Pars pars;
	pars.xsz = 100;
	pars.ysz = 4;
	pars.blx = C(-50.0);
	pars.bly = C(0.0);
	pars.dx = C(1.0);
	pars.dy = pars.dx;
	Pars *Parptr = &pars;

	Arrays arrays;
	Arrays *Arrptr = &arrays;

	Arrptr->DEM = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);
	Arrptr->H = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);
	dg2::allocate_fields(Parptr, Arrptr);

	dg2::initialise_field(Z, Parptr, Arrptr->DEM, Arrptr->DEM1x, Arrptr->DEM1y);
	dg2::initialise_h_from_eta(C(1.5), Parptr, Arrptr);

	dg2::write_dem("hump.dem", Parptr, Arrptr);
	dg2::write_startfile("hump.start", Parptr, Arrptr);

	dg2::deallocate_fields(Parptr, Arrptr);
	return EXIT_SUCCESS;
}
