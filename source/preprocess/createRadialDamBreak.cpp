#include "../lisflood.h"
#include "../utility.h"
#include "../swe/dg2/fields.h"
#include <cstdlib>
#include <cmath>

NUMERIC_TYPE H(NUMERIC_TYPE x, NUMERIC_TYPE y)
{
	NUMERIC_TYPE r = C(2.5), h_outside = C(0.5), h_inside = C(2.5);

	return (sqrt(x*x + y*y) <= r) ? h_inside : h_outside;
}

int main()
{
	Pars pars;
	pars.xsz = 200;
	pars.ysz = 200;
	pars.blx = C(-20.0);
	pars.tlx = C(20.0);
	pars.bly = C(-20.0);
	pars.tly = C(20.0);
	pars.dx = C(40.0) / pars.xsz;
	pars.dy = pars.dx;
	Pars *Parptr = &pars;

	Arrays arrays;
	Arrays *Arrptr = &arrays;

	Arrptr->DEM = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);
	Arrptr->H = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);
	dg2::allocate_fields(Parptr, Arrptr);

	auto Z = [](NUMERIC_TYPE x, NUMERIC_TYPE y) { return C(0.0); };

	dg2::initialise_field(Z, Parptr, Arrptr->DEM, Arrptr->DEM1x, Arrptr->DEM1y);
	dg2::initialise_field(H, Parptr, Arrptr->H, Arrptr->H1x, Arrptr->H1y);

	dg2::write_dem("radialDamBreak.dem", Parptr, Arrptr);
	dg2::write_startfile("radialDamBreak.start", Parptr, Arrptr);

	dg2::deallocate_fields(Parptr, Arrptr);
	return EXIT_SUCCESS;
}
