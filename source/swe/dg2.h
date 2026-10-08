#pragma once
#include "../lisflood.h"
#include "../rain/rain.h"

namespace dg2
{
	typedef struct
	{
		NUMERIC_TYPE *H;
		NUMERIC_TYPE *H1x;
		NUMERIC_TYPE *H1y;
		NUMERIC_TYPE *HU;
		NUMERIC_TYPE *HU1x;
		NUMERIC_TYPE *HU1y;
		NUMERIC_TYPE *HV;
		NUMERIC_TYPE *HV1x;
		NUMERIC_TYPE *HV1y;
	} FlowCoefficients;

	typedef struct
	{
		NUMERIC_TYPE H;
		NUMERIC_TYPE H1x;
		NUMERIC_TYPE H1y;
		NUMERIC_TYPE HU;
		NUMERIC_TYPE HU1x;
		NUMERIC_TYPE HU1y;
		NUMERIC_TYPE HV;
		NUMERIC_TYPE HV1x;
		NUMERIC_TYPE HV1y;
	} Increment;

    void zero_thin_depth_slopes
    (
        Pars*,
        Solver*,
        Arrays*
    );

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

	void update_point_sources
	(
		Pars *Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
		Arrays *Arrptr
	);

	void update_boundary_values
	(
		Pars *Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
		Arrays *Arrptr
	);

	void rk_stage1
	(
		Pars *Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
		Arrays *Arrptr
	);

	void rk_stage2
	(
		Pars *Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
		Arrays *Arrptr
	);

	void update_fluxes
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U
	);

	void update_fluxes_on_boundaries
	(
		Pars *Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
		Arrays *Arrptr,
		FlowCoefficients const& U
	);

	void set_boundary_values
	(
		Pars *Parptr,
		Solver *Solverptr,
		BoundCs *BCptr,
		const int bc_i,
		const int HU_sign,
		NUMERIC_TYPE Zstar,
		NUMERIC_TYPE& H_const,
		NUMERIC_TYPE HU_const,
		NUMERIC_TYPE HV_const,
		NUMERIC_TYPE& H_inside,
		NUMERIC_TYPE& HU_inside,
		NUMERIC_TYPE& HV_inside,
		NUMERIC_TYPE& H_outside,
		NUMERIC_TYPE& HU_outside,
		NUMERIC_TYPE& HV_outside
	);

	void update_intermediate_coefficients
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays* Arrptr,
		const int i,
		const int j,
		Increment const& U_inc
	);

	void update_intermediate_coefficient
	(
		Pars *Parptr,
		Solver *Solverptr,
		const int i,
		const int j,
		NUMERIC_TYPE *field_n,
		NUMERIC_TYPE *field_int,
		NUMERIC_TYPE increment
	);

	void update_final_coefficients
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays* Arrptr,
		const int i,
		const int j,
		Increment const& U_inc
	);

	void update_final_coefficient
	(
		Pars *Parptr,
		Solver *Solverptr,
		const int i,
		const int j,
		NUMERIC_TYPE *field_n,
		NUMERIC_TYPE *field_int,
		NUMERIC_TYPE increment
	);

	void set_initial_coefficients
	(
	 	Arrays *Arrptr,
		FlowCoefficients& U
	);

	void set_intermediate_coefficients
	(
	 	Arrays *Arrptr,
		FlowCoefficients& U_int
	);

	void L
	(
	 	Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j,
		FlowCoefficients const& U,
		Increment& U_inc
	);

	void L0
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j,
		FlowCoefficients const& U,
		Increment& U_inc
	);

	void L1x
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j,
		FlowCoefficients const& U,
		Increment& U_inc
	);

	void L1y
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j,
		FlowCoefficients const& U,
		Increment& U_inc
	);

	NUMERIC_TYPE bed_source_0x
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE bed_source_0y
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE bed_source_1x
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE bed_source_1y
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE maxH
	(
		Pars *Parptr,
		Arrays *Arrptr
	);

    bool thin_depth
    (
        Solver*,
        NUMERIC_TYPE H
    );

	NUMERIC_TYPE Tstep_from_cfl
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr
	);
	
	NUMERIC_TYPE min_dt
	(
		NUMERIC_TYPE dt,
		NUMERIC_TYPE H,
		NUMERIC_TYPE HU,
		NUMERIC_TYPE HV,
        Pars *Parptr,
		Solver *Solverptr
	);

    void zero_discharge
    (
        Pars*,
        Solver*,
        int i,
        int j,
        NUMERIC_TYPE* H,
        NUMERIC_TYPE* H1x,
        NUMERIC_TYPE* H1y,
        NUMERIC_TYPE* HU,
        NUMERIC_TYPE* HU1x,
        NUMERIC_TYPE* HU1y,
        NUMERIC_TYPE* HV,
        NUMERIC_TYPE* HV1x,
        NUMERIC_TYPE* HV1y
    );
}
