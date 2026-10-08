// LISFLOOD-FP Metal GPU back end: entry points called from lisflood.cpp
#pragma once
#include "../lisflood.h"

// true when the machine has a Metal GPU the FV1 solver can use
bool metal_is_available();

// run the FV1 solver on the Metal GPU
void run_metal_fv1(Fnames &filenames, Files *files, States &states, Pars &pars,
                   Solver &solver, BoundCs &boundCs, Stage &stage,
                   Arrays &arrays, int verbose);
