#include "../lisflood.h"
#include "../utility.h"
#include <cstdlib>
#include <cmath>

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

	NUMERIC_TYPE *DEM;
	NUMERIC_TYPE *H;

	DEM = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);
	H = memory_allocate_zero_numeric_legacy(pars.xsz * pars.ysz);

	NUMERIC_TYPE h_up = C(2.5), h_down = C(0.5);

#pragma omp parallel for
	for (int j=0; j<pars.ysz; j++)
	{
		for(int i=0; i<pars.xsz; i++)
		{
			NUMERIC_TYPE x = x_centre(&pars, i);
			NUMERIC_TYPE y = y_centre(&pars, j);
			H[j*pars.xsz + i] = x <= C(0.0) ? h_up : h_down;
		}
	}

	States states;
	states.alt_ascheader = OFF;
	states.call_gzip = OFF;

	write_ascfile("1DDamBreak", -1, ".dem", DEM, DEM, 0, &states, &pars);
	write_ascfile("1DDamBreak", -1, ".start", H, DEM, 0, &states, &pars);
	return EXIT_SUCCESS;
}
