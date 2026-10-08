#pragma once
#include "../lisflood.h"

void write_mass_stats
(
	Files *Fptr,
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr
);

void write_depth_samples
(
	Files *Fptr,
	Pars *Parptr,
	Solver *Solverptr,
	Stage *Stageptr,
	Arrays *Arrptr
);

void write_speed_samples
(
	Files *Fptr,
	Pars *Parptr,
	Solver *Solverptr,
	Stage *Stageptr,
	Arrays *Arrptr
);

void write_solution
(
	Fnames *Fnameptr,
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
);

void write_max_field
(
    Fnames*,
    States*,
    Pars*,
    Arrays*
);

void write_field
(
	NUMERIC_TYPE *field,
	const char *ascii_suffix,
	const char *binary_suffix,
	Fnames *Fnameptr,
	States *Statesptr,
	Pars *Parptr
);
