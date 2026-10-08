#include "dg2_output.h"
#include "../output.h"

void write_solution_slopes
(
	Fnames *Fnameptr,
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
    if (Statesptr->save_depth == ON)
	{
		write_field(Arrptr->H1x, ".wd1x", ".wd1xb", Fnameptr, Statesptr, Parptr);
		write_field(Arrptr->H1y, ".wd1y", ".wd1yb", Fnameptr, Statesptr, Parptr);
	}
	
    if (Statesptr->save_Qs == ON)
	{
		write_field(Arrptr->HU1x, ".Qx1x", ".Qx1xb", Fnameptr, Statesptr, Parptr);
		write_field(Arrptr->HU1y, ".Qx1y", ".Qx1yb", Fnameptr, Statesptr, Parptr);
		write_field(Arrptr->HV1x, ".Qy1x", ".Qy1xb", Fnameptr, Statesptr, Parptr);
		write_field(Arrptr->HV1y, ".Qy1y", ".Qy1yb", Fnameptr, Statesptr, Parptr);
	}
}
