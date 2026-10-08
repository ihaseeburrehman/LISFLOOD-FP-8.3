/*
* sgm_fast.h
*
*  Created on: 14 May 2014
*      Author: td14281
*/

#pragma once

#include "lisflood2.h"
#include "DataTypes.h"
#include "../rain/rain.h"

struct VolumeHeightUpdateInfo
{
	NUMERIC_TYPE evap_loss;
	NUMERIC_TYPE rain_total;
	NUMERIC_TYPE Qpoint_timestep_pos;
	NUMERIC_TYPE Qpoint_timestep_neg;
};

void SGC2_CalcA_public(int gr, NUMERIC_TYPE hflow, NUMERIC_TYPE bf, NUMERIC_TYPE *A, NUMERIC_TYPE *we, const SGCprams *SGCptr);

void SGC2_CalcLinksQ(SuperGridLinksList * Super_linksptr, NUMERIC_TYPE * volume_grid, const NUMERIC_TYPE * h_grid, const NUMERIC_TYPE delta_time, const NUMERIC_TYPE g, const NUMERIC_TYPE depth_thresh, const NUMERIC_TYPE max_Froude, WetDryRowBound * wet_dry_bounds);

NUMERIC_TYPE SGC2_CalculateVelocity_public(const int index, const int index_next,
	const NUMERIC_TYPE * Q_grid,
	const NUMERIC_TYPE * h_grid, const NUMERIC_TYPE * dem_grid, const NUMERIC_TYPE width);

void SGC2_UpdateLoadBalance(const int grid_rows, const int grid_cols_padded,
	const SubGridRowList * sub_grid_layout,
	WetDryRowBound* wet_dry_bounds);

void Fast_IterateLoop(const int grid_cols, const int grid_rows, const int grid_cols_padded,
	NUMERIC_TYPE *h_grid, NUMERIC_TYPE *volume_grid, 
	NUMERIC_TYPE *Qx_grid, NUMERIC_TYPE *Qy_grid, NUMERIC_TYPE *Qx_old_grid, NUMERIC_TYPE *Qy_old_grid,
	NUMERIC_TYPE *maxH_grid, NUMERIC_TYPE *maxHtm_grid, NUMERIC_TYPE *initHtm_grid, NUMERIC_TYPE *totalHtm_grid,
	NUMERIC_TYPE *maxVc_grid, NUMERIC_TYPE *maxVc_height_grid, NUMERIC_TYPE *maxHazard_grid,
	NUMERIC_TYPE *Vx_grid, NUMERIC_TYPE *Vy_grid, NUMERIC_TYPE *Vx_max_grid, NUMERIC_TYPE *Vy_max_grid,
	const NUMERIC_TYPE *dem_grid,
	const NUMERIC_TYPE *g_friction_sq_x_grid, const NUMERIC_TYPE *g_friction_sq_y_grid,
	const NUMERIC_TYPE *friction_x_grid, const NUMERIC_TYPE *friction_y_grid,
	const NUMERIC_TYPE *dx_col, const NUMERIC_TYPE *dy_col, const NUMERIC_TYPE *cell_area_col,
	const NUMERIC_TYPE *Fp_xwidth, const NUMERIC_TYPE *Fp_ywidth,

	const SubGridRowList * sub_grid_layout_rows,
	SubGridState * sub_grid_state_rows,
	const SubGridRowList * sub_grid_layout_blocks,
	SubGridState * sub_grid_state_blocks,

	const NUMERIC_TYPE * SGC_BankFullHeight_grid,

	TimeSeries * evap_time_series,
	NetCDFVariable *evap_grid,	      
	TimeSeries * rain_time_series,
	NUMERIC_TYPE *rain_grid,
	const NUMERIC_TYPE *dist_infil_grid,

	WetDryRowBound* wet_dry_bounds,
	PointSourceRowList * ps_layout, BoundaryCondition * boundary_cond,
	WeirLayout * weirs_weirs, WeirLayout * weirs_bridges,
	RouteDynamicList * route_dynamic_list,
	const NUMERIC_TYPE *route_V_ratio_per_sec_qx, const NUMERIC_TYPE * route_V_ratio_per_sec_qy,

	SuperGridLinksList * Super_linksptr,
	Fnames *Fnameptr,
	Files *Fptr,
	Stage *Locptr,
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	DamData *Damptr,
	SGCprams * SGCptr,
	NUMERIC_TYPE ** tmp_thread_data,
#ifdef RESULT_CHECK
	Arrays * Arrptr, // only for compare results
	BoundCs * BCptr, // only for compare results
	ChannelSegmentType *ChannelSegments,
	vector<ChannelSegmentType> *ChannelSegmentsVecPtr,
#endif

	const int verbose);
