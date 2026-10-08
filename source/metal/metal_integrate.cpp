// LISFLOOD-FP Metal GPU back end: entry points called from lisflood.cpp
#include "metal_integrate.h"
#include "fv1/metal_fv1_simulate.h"

bool metal_is_available() { return lis::metal::fv1::Simulator::is_available(); }

void run_metal_fv1(Fnames &filenames, Files *files, States &states, Pars &pars,
                   Solver &solver, BoundCs &boundCs, Stage &stage,
                   Arrays &arrays, int verbose) {
  if (verbose)
    printf("\nUsing the Metal GPU (Apple Silicon) for the FV1 solver\n\n");
  lis::metal::fv1::simulate(filenames, files, states, pars, solver, boundCs,
                            stage, arrays, verbose);
}
