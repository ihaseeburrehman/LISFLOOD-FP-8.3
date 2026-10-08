#include "../lisflood.h"
#include "../utility.h"
#include "../swe/dg2/fields.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <initializer_list>

NUMERIC_TYPE Z(NUMERIC_TYPE x, NUMERIC_TYPE y)
{
	NUMERIC_TYPE a = C(1.0) - C(0.2)*sqrt(pow(x-C(20.0),2) + pow(y-C(15.0),2));
	NUMERIC_TYPE b = C(2.0) - C(0.5)*sqrt(pow(x-C(40.0),2) + pow(y-C(15.0),2));
	NUMERIC_TYPE c = C(1.78) - C(0.3)*sqrt(pow(x-C(60.0),2) + pow(y-C(15.0),2));

	return std::max({C(0.0), a, b, c});
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
	dg2::initialise_h_from_eta(C(1.78), Parptr, Arrptr);

	dg2::write_dem("lakeAtRestCones.dem", Parptr, Arrptr);
	dg2::write_startfile("lakeAtRestCones.start", Parptr, Arrptr);

	dg2::deallocate_fields(Parptr, Arrptr);
	return EXIT_SUCCESS;
}
