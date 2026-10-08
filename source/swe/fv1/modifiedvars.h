#pragma once
#include "../../lisflood.h"

namespace fv1
{
	void initialise_Zstar
	(
		Pars *Parptr,
		Arrays *Arrptr
	);

	void update_Hstar
	(
		Pars *Parptr,
		Arrays *Arrptr
	);

	NUMERIC_TYPE HUstar_neg_x
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HUstar_pos_x
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HUstar_neg_y
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HUstar_pos_y
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HVstar_neg_x
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HVstar_pos_x
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HVstar_neg_y
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE HVstar_pos_y
	(
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE discharge_star
	(
		NUMERIC_TYPE *Hstar,
		NUMERIC_TYPE *discharge_component,
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger_neg_x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger_pos_x
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger_neg_y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger_pos_y
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE Zdagger
	(
		Pars *Parptr,
		Arrays *Arrptr,
		NUMERIC_TYPE *Zstar,
		const int Zstar_i,
		const int Zstar_j,
		const int ETA_i,
		const int ETA_j
	);

	NUMERIC_TYPE eta
	(
		Pars *Parptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);

	NUMERIC_TYPE speed
	(
		NUMERIC_TYPE *discharge_component,
		Pars *Parptr,
		Solver *Solverptr,
		Arrays *Arrptr,
		const int i,
		const int j
	);
}
