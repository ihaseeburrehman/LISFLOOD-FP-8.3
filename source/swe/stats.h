#pragma once
#include "../lisflood.h"

void update_mass_stats
(
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr
);

NUMERIC_TYPE flood_area
(
	Pars *Parptr,
	Arrays *Arrptr
);

void zero_flux_stats(BoundCs *BCptr);

void accumulate_point_flux_stats
(
	BoundCs *BCptr,
	NUMERIC_TYPE Q
);

void accumulate_boundary_flux_stats
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr
);

void update_velocity
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
);

void update_max_field(Pars*, Arrays*);
