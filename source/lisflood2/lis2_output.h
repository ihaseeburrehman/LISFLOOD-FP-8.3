#pragma once

#include "../lisflood.h"
#include "lisflood2.h"
#include "DataTypes.h"

void write_regular_output(const char * resrootname, const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE curr_time,
	NUMERIC_TYPE*  tmp_grid1, NUMERIC_TYPE*  tmp_grid2, NUMERIC_TYPE*  tmp_grid3,
	const NUMERIC_TYPE* h_grid, const NUMERIC_TYPE* dem_grid,
	const NUMERIC_TYPE* Qx_grid, const NUMERIC_TYPE* Qy_grid,
	const NUMERIC_TYPE* Qx_old_grid, const NUMERIC_TYPE* Qy_old_grid,
	const NUMERIC_TYPE* Vx_grid, const NUMERIC_TYPE* Vy_grid,
	const NUMERIC_TYPE* SGC_BankFullHeight_grid, const NUMERIC_TYPE* maxH_grid,
	const SubGridRowList * sub_grid_layout, const SubGridState * sub_grid_state,

	const NUMERIC_TYPE* dx_col, const NUMERIC_TYPE* dy_col,
	Solver *Solverptr, States *Statesptr, Pars *Parptr, BoundaryCondition *boundary_cond, SGCprams *SGCptr,
	OutputParams * output_params,
	const int save_number,
	const int save_depth,
	const int save_elev,
	const int save_Qs,
	const int save_Velocity,
	int save_sub_grid_Velocity);

void WriteOutput(Fnames *Fnameptr, const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh,
	NUMERIC_TYPE*  tmp_grid,
	const NUMERIC_TYPE * initHtm_grid, const NUMERIC_TYPE * totalHtm_grid, const NUMERIC_TYPE * maxH_grid, const NUMERIC_TYPE * maxHtm_grid,
	const NUMERIC_TYPE * maxVc_grid, const NUMERIC_TYPE * maxVc_height_grid, const NUMERIC_TYPE * maxHazard_grid,
	const NUMERIC_TYPE * Vx_max_grid, const NUMERIC_TYPE * Vy_max_grid,
	const NUMERIC_TYPE * dem_grid,
	const NUMERIC_TYPE * SGC_BankFullHeight_grid,
	States *Statesptr, Pars *Parptr, OutputParams* output_params);