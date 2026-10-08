#pragma once
#include "../dg2.h"
#include "../../lisflood.h"

namespace dg2
{
	void initialise_Zstar
	(
		Pars *Parptr,
		Arrays *Arrptr
	);

	void update_Hstar
	(
		Pars *Parptr,
		Solver* Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U
	);

	void update_HUstar
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U
	);

	void update_HVstar
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		FlowCoefficients const& U
	);

	void update_discharge_star
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		NUMERIC_TYPE* discharge_const,
		NUMERIC_TYPE* discharge1x,
		NUMERIC_TYPE* discharge1y,
		FlowCoefficients const& U,
		NUMERIC_TYPE* discharge_star_neg_x,
		NUMERIC_TYPE* discharge_star_pos_x,
		NUMERIC_TYPE* discharge_star_neg_y,
		NUMERIC_TYPE* discharge_star_pos_y
	);

	NUMERIC_TYPE discharge_star
	(
		Solver *Solverptr,
		NUMERIC_TYPE H,
		NUMERIC_TYPE discharge,
		NUMERIC_TYPE Hstar
	);

	NUMERIC_TYPE Hstar0x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE Hstar0y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE Hstar1x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE Hstar1y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HUstar0x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HUstar0y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HUstar1x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HUstar1y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HVstar0x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HVstar0y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HVstar1x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HVstar1y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger1x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger1y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger_neg_x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger_pos_x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger_neg_y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger_pos_y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE eta_neg_x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE eta_pos_x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE eta_neg_y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE eta_pos_y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE eta
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE eta1x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE eta1y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		FlowCoefficients const& U,
		const int i,
		const int j
	);

	NUMERIC_TYPE limit_neg
	(
		NUMERIC_TYPE *constant,
		NUMERIC_TYPE *slope,
		Pars *Parptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE limit_pos
	(
		NUMERIC_TYPE *constant,
		NUMERIC_TYPE *slope,
		Pars *Parptr,
		const int i,
		const int j
	);
	
	NUMERIC_TYPE gauss_lower
	(
		NUMERIC_TYPE *constant,
		NUMERIC_TYPE *slope,
		Pars *Parptr,
		const int i,
		const int j
	);
	
	NUMERIC_TYPE gauss_upper
	(
		NUMERIC_TYPE *constant,
		NUMERIC_TYPE *slope,
		Pars *Parptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE limit_neg(NUMERIC_TYPE constant, NUMERIC_TYPE slope);
	NUMERIC_TYPE limit_pos(NUMERIC_TYPE constant, NUMERIC_TYPE slope);
	NUMERIC_TYPE gauss_lower(NUMERIC_TYPE constant, NUMERIC_TYPE slope);
	NUMERIC_TYPE gauss_upper(NUMERIC_TYPE constant, NUMERIC_TYPE slope);
}

