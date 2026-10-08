#pragma once
#include "../lisflood.h"
#include "../rain/rain.h"

namespace fv1
{
	void solve
	(
		Fnames *Fnameptr,
		Files *Fptr,
		States *Statesptr,
		Pars *Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
		Stage *Stageptr,
		Arrays *Arrptr,
		const int verbose
	);
    
    void drain_nodata_water
    (
        Pars* Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
        Arrays *Arrptr
    );

    void update_point_sources
    (
        Pars *Parptr,
        Solver *Solverptr,
        BoundCs *BCptr,
        Arrays *Arrptr
    );

	void apply_friction
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr
	);

	void update_fluxes
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr
	);

	void update_fluxes_on_boundaries
	(
		Pars *Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
		Arrays *Arrptr
	);

	void update_flow_variables
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr
	);

	void set_boundary_values
	(
		Pars *Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
		Arrays *Arrptr,
		const int z_i,
		const int z_j,
		int bc_i,
		int HU_sign,
		NUMERIC_TYPE& H_inside,
		NUMERIC_TYPE& HU_inside,
		NUMERIC_TYPE& HV_inside,
		NUMERIC_TYPE& H_outside,
		NUMERIC_TYPE& HU_outside,
		NUMERIC_TYPE& HV_outside
	);

	NUMERIC_TYPE Tstep_from_cfl
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr
	);

	NUMERIC_TYPE bed_source_x
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE bed_source_y
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);
}
