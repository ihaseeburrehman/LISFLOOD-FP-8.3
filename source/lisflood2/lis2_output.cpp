#include "lis2_output.h"
#include "file_tool.h"
#include "sgm_fast.h"

/*
grid_cols_padded required padding of destination grid
*/
void add_boundary_sub_flow(const int grid_cols, const int grid_rows, const int grid_cols_padded,
	NUMERIC_TYPE*  tmp_grid,
	BoundaryCondition * boundary_cond,
	const bool is_qx // (otherwise qy)
	)
{
	if (!is_qx)
	{
		// North Boundary
		for (int i = 0; i < grid_cols; i++)
		{
			int BCi = i;
			NUMERIC_TYPE q = boundary_cond->bc_info.Q_SG_old[BCi]; //note SGM_Qy_grid will always be 0 in case of _linear sub grid
			tmp_grid[i] += q;

		}
		// South boundary
		for (int i = 0; i < grid_cols; i++)
		{
			int BCi = i + (grid_cols + grid_rows);
			int x = grid_cols - i - 1;
			int q_index = x + grid_rows * (grid_cols_padded);
			NUMERIC_TYPE q = boundary_cond->bc_info.Q_SG_old[BCi];
			tmp_grid[q_index] += q;
		}
	}
	else
	{
		//	West boundary
		for (int j = 0; j < grid_rows; j++)
		{
			int BCi = j + (2 * grid_cols + grid_rows);
			int y = grid_rows - j - 1;
			int q_index = y * (grid_cols_padded);
			NUMERIC_TYPE q = boundary_cond->bc_info.Q_SG_old[BCi];
			tmp_grid[q_index] += q;
		}

		// East boundary
		for (int j = 0; j < grid_rows; j++)
		{
			int BCi = grid_cols + j;
			int q_index = j * grid_cols_padded + grid_cols;
			NUMERIC_TYPE q = boundary_cond->bc_info.Q_SG_old[BCi];
			tmp_grid[q_index] += q;
		}
	}
}

/*
takes the values in the 'flow_list' array and puts into the correct grid location
flow_list could be:
* sub_grid_state->sg_flow_Q or
* sub_grid_state->sg_velocity
* sub_grid_state->sg_flow_ChannelRatio

grid_cols_padded required padding of destination grid
*/
void get_sub_grid_values(const int grid_cols, const int grid_rows, const int grid_cols_padded,
	NUMERIC_TYPE*  dest_grid,
	const NUMERIC_TYPE * source_flow_list,
	const SubGridRowList * sub_grid_layout,
	const bool is_qx // (otherwise qy)
	)
{
	SubGridCellInfo flow_pair = sub_grid_layout->flow_info.flow_pair;

	const int row_cols_padded = sub_grid_layout->row_cols_padded;

	for (int j = 0; j < grid_rows; j++)
	{
		int sg_row_start_index = j * row_cols_padded;
		int sg_row_pair_start_index = j * 2 * row_cols_padded;

		const int flow_end = sub_grid_layout->flow_row_count[j];
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
		__assume(sg_row_start_index % GRID_ALIGN_WIDTH == 0);
		__assume(sg_row_pair_start_index % GRID_ALIGN_WIDTH == 0);
#endif
#pragma ivdep
		for (int flow_i = 0; flow_i < flow_end; flow_i++)
		{
			int flow_index = sg_row_start_index + flow_i;
			int flow_pair_index = sg_row_pair_start_index + 2 * flow_i;
			int flow_pair_index_next = flow_pair_index + 1;



			//int grid_index = flow_pair.sg_cell_grid_index_lookup[flow_pair_index];
			//int grid_index_next = flow_pair.sg_cell_grid_index_lookup[flow_pair_index_next];
			int grid_index = flow_pair.sg_cell_x[flow_pair_index] + flow_pair.sg_cell_y[flow_pair_index] * grid_cols_padded;
			int grid_index_next = flow_pair.sg_cell_x[flow_pair_index_next] + flow_pair.sg_cell_y[flow_pair_index_next] * grid_cols_padded;

			if (is_qx)
			{
				if ((grid_index + 1) == grid_index_next)
					dest_grid[grid_index_next] = source_flow_list[flow_index];
			}
			else
			{
				if ((grid_index + 1) != grid_index_next)
					dest_grid[grid_index_next] = source_flow_list[flow_index];
			}
		}
	}

}


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
	const int save_sub_grid_Velocity)
{
	if (Statesptr->maxdepthonly == ON)
	{
	}
	else
	{
		
			int save_sub_grid;
		// if no sub-grid cells, then row_cols_padded will be zero -> disable sub grid writing
		if (sub_grid_layout->row_cols_padded > 0)
		{
			save_sub_grid = ON;
		}
		else
		{
			save_sub_grid = OFF;
		}

		if (output_params->netcdf_out == ON && output_params->netcdf_state.init_done == OFF)
		{
			char netcdf_filename[512];
			sprintf(netcdf_filename, "%s.nc", resrootname);

			StartNetCDF(grid_cols, grid_rows, netcdf_filename, &output_params->netcdf_state, save_depth, save_elev, save_Qs, save_Velocity, save_sub_grid_Velocity, save_sub_grid);
		}
		write_time(save_number, curr_time, output_params);

		// write out water depth
		if (save_depth == ON)
		{
			// update h to be 'depth' above sub-grid-channel
			for (int j = 0; j < grid_rows; j++)
			{
				for (int i = 0; i < grid_cols; i++)
				{
					int index_padded = i + j * grid_cols_padded;
					int index = i + j * grid_cols;
					NUMERIC_TYPE temp = h_grid[index_padded] + SGC_BankFullHeight_grid[index_padded];
					// depth in channel or above flood plain = h + BankFullHeight (should not be negative)
					// note previous version just dumped the wd with no depth_thresh truncation
					if (temp <= depth_thresh)
						temp = C(0.0);
					tmp_grid1[index] = temp;
				}
			}
			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_depth, ".wd",
				tmp_grid1, grid_cols, grid_rows, Parptr->blx, Parptr->bly, Parptr->dx, output_params);

			// export floodplain depth for SGC model only
			// h is height above the flood-plain - set to zero if below the flood plain
			for (int j = 0; j < grid_rows; j++)
			{
				for (int i = 0; i < grid_cols; i++)
				{
					int index_padded = i + j * grid_cols_padded;
					int index = i + j * grid_cols;
					// calculate depth on floodplain
					NUMERIC_TYPE temp = h_grid[index_padded];
					// if water level is below floodplain set depth to zero
					if (temp <= depth_thresh)
						temp = C(0.0);
					tmp_grid1[index] = temp;
				}
			}
			write_grid(resrootname, save_number,
				NETCDF_IGNORE, ".wdfp",
				tmp_grid1, grid_cols, grid_rows, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
		}

		// write out maximum water depth at each save interval and reset the values
		if (save_depth == ON && Statesptr->saveint_max == ON)
		{
			// update h to be 'depth' above sub-grid-channel
			for (int j = 0; j < grid_rows; j++)
			{
				for (int i = 0; i < grid_cols; i++)
				{
					int index_padded = i + j * grid_cols_padded;
					int index = i + j * grid_cols;
					NUMERIC_TYPE temp = maxH_grid[index_padded] + SGC_BankFullHeight_grid[index_padded];
					// depth in channel or above flood plain = h + BankFullHeight (should not be negative)
					// note previous version just dumped the wd with no depth_thresh truncation
					if (temp <= depth_thresh)
						temp = C(0.0);
					tmp_grid1[index] = temp;
				}
			}
			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_depth, ".wd_max",
				tmp_grid1, grid_cols, grid_rows, Parptr->blx, Parptr->bly, Parptr->dx, output_params);

			// export floodplain depth for SGC model only
			// h is height above the flood-plain - set to zero if below the flood plain
			for (int j = 0; j < grid_rows; j++)
			{
				for (int i = 0; i < grid_cols; i++)
				{
					int index_padded = i + j * grid_cols_padded;
					int index = i + j * grid_cols;
					// calculate depth on floodplain
					NUMERIC_TYPE temp = maxH_grid[index_padded];
					// if water level is below floodplain set depth to zero
					if (temp <= depth_thresh)
						temp = C(0.0);
					tmp_grid1[index] = temp;
				}
			}
			write_grid(resrootname, save_number,
				NETCDF_IGNORE, ".wdfp_max",
				tmp_grid1, grid_cols, grid_rows, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
		}

		// write out elevation data
		if (save_elev == ON)
		{
			for (int j = 0; j < grid_rows; j++)
			{
				for (int i = 0; i < grid_cols; i++)
				{
					int index_padded = i + j * grid_cols_padded;
					int index = i + j * grid_cols;
					// calculate depth on floodplain
					NUMERIC_TYPE temp = h_grid[index_padded];
					// if water level is below floodplain or sub-grid-channel set depth to NULLVAL
					if (temp + SGC_BankFullHeight_grid[index_padded] <= depth_thresh)
					{
						temp = NULLVAL;
					}
					else
					{
						temp += dem_grid[index_padded];
					}
					tmp_grid1[index] = temp;
				}
			}

			// if using SGC the elevation neds to be from the channel bed
			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_elevation, ".elev",
				tmp_grid1, grid_cols, grid_rows, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
		}

		const int qx_grid_cols = (grid_cols + 1); // columns without padding
		const int qy_grid_cols = (grid_cols); // columns without padding

		// write out x&y "flow fluxes" if required (note flag=1 for Qx and flag=2 for Qy)
		if (save_Qs == ON)
		{
			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_qx, ".Qx",
				Qx_grid, tmp_grid1, qx_grid_cols, grid_rows, grid_cols_padded, Parptr->blx - (Parptr->dx / C(2.0)), Parptr->bly, Parptr->dx, output_params);

			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_qy, ".Qy",
				Qy_grid, tmp_grid1, qy_grid_cols, grid_rows + 1, grid_cols_padded, Parptr->blx, Parptr->bly - (Parptr->dy / C(2.0)), Parptr->dx, output_params);

			if (save_sub_grid == ON)
			{
				memset(tmp_grid1, 0, sizeof(NUMERIC_TYPE) * (grid_rows + 1) * grid_cols_padded);
				// additional sub grid output of flow width
				for (int j = 0; j < grid_rows; j++)
				{
					int sg_list_row_start = j * sub_grid_layout->row_cols_padded;
#pragma ivdep
					for (int cell_i = 0; cell_i < sub_grid_layout->cell_row_count[j]; cell_i++)
					{
						int cell_index = cell_i + sg_list_row_start;
						int grid_index_padded = sub_grid_layout->cell_info.sg_cell_grid_index_lookup[cell_index];
						int grid_index = sub_grid_layout->cell_info.sg_cell_x[cell_index] + sub_grid_layout->cell_info.sg_cell_y[cell_index] * grid_cols;

						NUMERIC_TYPE A0, w0;
						A0 = C(0.0);
						NUMERIC_TYPE bank_full_height = sub_grid_layout->cell_info.sg_cell_SGC_BankFullHeight[cell_index];
						NUMERIC_TYPE h = h_grid[grid_index_padded];
						h += bank_full_height;
						if (h > C(0.0))
						{
							w0 = sub_grid_layout->cell_info.sg_cell_SGC_width[cell_index];
							SGC2_CalcA_public(sub_grid_layout->cell_info.sg_cell_SGC_group[cell_index], h, bank_full_height, &A0, &w0, SGCptr);
						}
						else
						{
							w0 = C(0.0);
						}
						tmp_grid1[grid_index] = w0;
					}
				}

				write_grid(resrootname, save_number,
					NETCDF_IGNORE, ".Fwidth",
					tmp_grid1, grid_cols, grid_rows, Parptr->blx, Parptr->bly, Parptr->dx, output_params);

				memset(tmp_grid1, 0, sizeof(NUMERIC_TYPE) * (grid_rows + 1) * grid_cols_padded);
				get_sub_grid_values(grid_cols, grid_rows, qx_grid_cols, tmp_grid1, sub_grid_state->sg_flow_Q, sub_grid_layout, true);
				add_boundary_sub_flow(grid_cols, grid_rows, qx_grid_cols, tmp_grid1, boundary_cond, true);
				write_grid(resrootname, save_number,
					output_params->netcdf_state.varid_qcx, ".Qcx",
					tmp_grid1, qx_grid_cols, grid_rows, Parptr->blx - (Parptr->dx / C(2.0)), Parptr->bly, Parptr->dx, output_params);

				memset(tmp_grid1, 0, sizeof(NUMERIC_TYPE) * (grid_rows + 1) * grid_cols_padded);
				get_sub_grid_values(grid_cols, grid_rows, qy_grid_cols, tmp_grid1, sub_grid_state->sg_flow_Q, sub_grid_layout, false);
				add_boundary_sub_flow(grid_cols, grid_rows, qy_grid_cols, tmp_grid1, boundary_cond, false);
				write_grid(resrootname, save_number,
					output_params->netcdf_state.varid_qcy, ".Qcy",
					tmp_grid1, qy_grid_cols, grid_rows + 1, Parptr->blx, Parptr->bly - (Parptr->dy / C(2.0)), Parptr->dx, output_params);
			}
		}

		// write out x&y adaptive timesteps if required
		if (Statesptr->voutput == ON)
		{
			memset(tmp_grid1, 0, sizeof(NUMERIC_TYPE) * (grid_rows + 1) * grid_cols_padded);

			// copy removing padding
			for (int j = 0; j < grid_rows; j++)
			{
				memcpy(tmp_grid1 + (j * qx_grid_cols), Vx_grid + (j * grid_cols_padded), sizeof(NUMERIC_TYPE) * qx_grid_cols);
			}
			/// cells with large sub grid do not have a flood plain flow,
			/// -> so reset the velocity to zero
			for (int j = 0; j < grid_rows; j++)
			{
				int sg_list_row_start = j * sub_grid_layout->row_cols_padded;
#pragma ivdep
				for (int cell_i = 0; cell_i < sub_grid_layout->cell_row_count[j]; cell_i++)
				{
					int cell_index = cell_i + sg_list_row_start;
					if (sub_grid_layout->cell_info.sg_cell_SGC_is_large[cell_index])
					{
						int grid_index = sub_grid_layout->cell_info.sg_cell_x[cell_index] +
							sub_grid_layout->cell_info.sg_cell_y[cell_index] * qx_grid_cols;
						tmp_grid1[grid_index] = C(0.0);
					}
				}
			}
			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_Vx, ".Vx",
				tmp_grid1, qx_grid_cols, grid_rows, Parptr->blx - (Parptr->dx / C(2.0)), Parptr->bly, Parptr->dx, output_params);

			memset(tmp_grid1, 0, sizeof(NUMERIC_TYPE) * (grid_rows + 1) * grid_cols_padded);
			for (int j = 0; j < grid_rows - 1; j++)
			{
				memcpy(tmp_grid1 + (j * qy_grid_cols), Vy_grid + j * grid_cols_padded, sizeof(NUMERIC_TYPE) * qy_grid_cols);
			}
			/// cells with large sub grid do not have a flood plain flow, so reset the velocity to zero
			for (int j = 0; j < grid_rows; j++)
			{
				int sg_list_row_start = j * sub_grid_layout->row_cols_padded;
#pragma ivdep
				for (int cell_i = 0; cell_i < sub_grid_layout->cell_row_count[j]; cell_i++)
				{
					int cell_index = cell_i + sg_list_row_start;
					if (sub_grid_layout->cell_info.sg_cell_SGC_is_large[cell_index])
					{
						int grid_index = sub_grid_layout->cell_info.sg_cell_x[cell_index] +
							sub_grid_layout->cell_info.sg_cell_y[cell_index] * qy_grid_cols;
						tmp_grid1[grid_index] = C(0.0);
					}
				}
			}
			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_Vy, ".Vy",
				tmp_grid1, qy_grid_cols, grid_rows + 1, Parptr->blx, Parptr->bly - (Parptr->dy / C(2.0)), Parptr->dx, output_params);
		}
		// output SGC channel velocity
		if (save_sub_grid == ON && Statesptr->SGCvoutput == ON)
		{
			int qx_grid_cols = (grid_cols + 1);
			memset(tmp_grid1, 0, sizeof(NUMERIC_TYPE) * (grid_rows + 1) * grid_cols_padded);
			get_sub_grid_values(grid_cols, grid_rows, qx_grid_cols, tmp_grid1, sub_grid_state->sg_velocity, sub_grid_layout, true);

			int qy_grid_cols = (grid_cols);
			memset(tmp_grid2, 0, sizeof(NUMERIC_TYPE) * (grid_rows + 1) * grid_cols_padded);
			get_sub_grid_values(grid_cols, grid_rows, qy_grid_cols, tmp_grid1, sub_grid_state->sg_velocity, sub_grid_layout, false);

			memset(tmp_grid3, 0, sizeof(NUMERIC_TYPE) * (grid_rows + 1) * grid_cols_padded);
			//// Get largest velocity from cell faces and place at cell centre.
			for (int j = 0; j < grid_rows; j++)
			{
				int sg_list_row_start = j * sub_grid_layout->row_cols_padded;
#pragma ivdep
				for (int cell_i = 0; cell_i < sub_grid_layout->cell_row_count[j]; cell_i++)
				{
					int cell_index = cell_i + sg_list_row_start;
					if (sub_grid_layout->cell_info.sg_cell_SGC_is_large[cell_index])
					{
						int grid_index_qx = sub_grid_layout->cell_info.sg_cell_x[cell_index] +
							sub_grid_layout->cell_info.sg_cell_y[cell_index] * qx_grid_cols;
						int grid_index_qy = sub_grid_layout->cell_info.sg_cell_x[cell_index] +
							sub_grid_layout->cell_info.sg_cell_y[cell_index] * qy_grid_cols;
						int grid_index_out = sub_grid_layout->cell_info.sg_cell_x[cell_index] +
							sub_grid_layout->cell_info.sg_cell_y[cell_index] * grid_cols;

						tmp_grid3[grid_index_out] =
							getmax(
								getmax(FABS(tmp_grid1[grid_index_qx]), FABS(tmp_grid1[grid_index_qx + 1])),
								getmax(FABS(tmp_grid2[grid_index_qy]), FABS(tmp_grid2[grid_index_qy + qy_grid_cols]))
								);
					}
				}
			}
			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_sgc_Vx, ".SGCVx",
				tmp_grid1, qx_grid_cols, grid_rows + 1, Parptr->blx - (Parptr->dx / C(2.0)), Parptr->bly, Parptr->dx, output_params);
			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_sgc_Vy, ".SGCVy",
				tmp_grid2, qy_grid_cols, grid_rows + 1, Parptr->blx, Parptr->bly - (Parptr->dy / C(2.0)), Parptr->dx, output_params);
			write_grid(resrootname, save_number,
				output_params->netcdf_state.varid_sgc_Vc, ".SGCVc",
				tmp_grid3, grid_cols, grid_rows + 1, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
		}
		// write out x&y Q limits if required
		if (Statesptr->save_QLs == ON)
		{
			//todo
			//if (Statesptr->binary_out == ON) // output binary of ascii rasters
			//{
			//	write_binrasterfile(resrootname, save_number, ".QLxb", Arrptr->LimQx, Arrptr->DEM, 1, Statesptr, Parptr);
			//	write_binrasterfile(resrootname, save_number, ".QLyb", Arrptr->LimQy, Arrptr->DEM, 2, Statesptr, Parptr);
			//}
			//else
			//{
			//	write_ascfile(resrootname, save_number, ".QLx", Arrptr->LimQx, Arrptr->DEM, 1, Statesptr, Parptr);
			//	write_ascfile(resrootname, save_number, ".QLy", Arrptr->LimQy, Arrptr->DEM, 2, Statesptr, Parptr);
			//}
		}
	}

	return;
}

void WriteOutput(Fnames *Fnameptr, const int grid_cols, const int grid_rows, const int grid_cols_padded,
	const NUMERIC_TYPE depth_thresh,
	NUMERIC_TYPE*  tmp_grid,
	const NUMERIC_TYPE * initHtm_grid, const NUMERIC_TYPE * totalHtm_grid, const NUMERIC_TYPE * maxH_grid, const NUMERIC_TYPE * maxHtm_grid,
	const NUMERIC_TYPE * maxVc_grid, const NUMERIC_TYPE * maxVc_height_grid, const NUMERIC_TYPE * maxHazard_grid,
	const NUMERIC_TYPE * Vx_max_grid, const NUMERIC_TYPE * Vy_max_grid,
	const NUMERIC_TYPE * dem_grid,
	const NUMERIC_TYPE * SGC_BankFullHeight_grid,
	States *Statesptr, Pars *Parptr, OutputParams* output_params)
{
	if (Statesptr->maxdepthonly== ON)
	{ 	// Write maximum depth
	write_grid(Fnameptr->resrootname, -1,
		NETCDF_IGNORE,
		".max", maxH_grid, tmp_grid, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, output_params);

	for (int j = 0; j < grid_rows; j++)
	{
		for (int i = 0; i < grid_cols; i++)
		{
			int index_padded = i + j * grid_cols_padded;
			int index = i + j * grid_cols;
			NUMERIC_TYPE temp;
			if (maxH_grid[index_padded] <= depth_thresh)
				temp = NULLVAL;
			else
			{
				temp = maxH_grid[index_padded] + dem_grid[index_padded] - SGC_BankFullHeight_grid[index_padded];
			}
			tmp_grid[index] = temp;
		}
	}
	}
	else
	{
		// Write time of initial flood inundation
		write_grid(Fnameptr->resrootname, -1,
			NETCDF_IGNORE, ".inittm",
			initHtm_grid, tmp_grid, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
		// Write total inundation time
		write_grid(Fnameptr->resrootname, -1,
			NETCDF_IGNORE, ".totaltm",
			totalHtm_grid, tmp_grid, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
		// Write maximum depth
		write_grid(Fnameptr->resrootname, -1,
			NETCDF_IGNORE,
			".max", maxH_grid, tmp_grid, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, output_params);

		for (int j = 0; j < grid_rows; j++)
		{
			for (int i = 0; i < grid_cols; i++)
			{
				int index_padded = i + j * grid_cols_padded;
				int index = i + j * grid_cols;
				NUMERIC_TYPE temp;
				if (maxH_grid[index_padded] <= depth_thresh)
					temp = NULLVAL;
				else
				{
					temp = maxH_grid[index_padded] + dem_grid[index_padded] - SGC_BankFullHeight_grid[index_padded];
				}
				tmp_grid[index] = temp;
			}
		}
		//// Write maximum elevation
		write_grid(Fnameptr->resrootname, -1,
			NETCDF_IGNORE, ".mxe",
			tmp_grid, grid_cols, grid_rows, Parptr->blx, Parptr->bly, Parptr->dx, output_params);

		// Write time of maximum depth
		write_grid(Fnameptr->resrootname, -1,
			NETCDF_IGNORE, ".maxtm",
			maxHtm_grid, tmp_grid, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
		if (Statesptr->voutput == ON)
		{
			write_grid(Fnameptr->resrootname, -1,
				NETCDF_IGNORE, ".maxVx",
				Vx_max_grid, tmp_grid, grid_cols + 1, grid_rows, grid_cols_padded, Parptr->blx - (Parptr->dx / C(2.0)), Parptr->bly, Parptr->dx, output_params);
			write_grid(Fnameptr->resrootname, -1,
				NETCDF_IGNORE, ".maxVy",
				Vy_max_grid, tmp_grid, grid_cols, grid_rows + 1, grid_cols_padded, Parptr->blx, Parptr->bly - (Parptr->dy / C(2.0)), Parptr->dx, output_params);
		}
		if (Statesptr->hazard == ON)
		{
			// Write maximum V Vd and Hazard
			write_grid(Fnameptr->resrootname, -1,
				NETCDF_IGNORE, ".maxVc",
				maxVc_grid, tmp_grid, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
			write_grid(Fnameptr->resrootname, -1,
				NETCDF_IGNORE, ".maxVcd",
				maxVc_height_grid, tmp_grid, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
			write_grid(Fnameptr->resrootname, -1,
				NETCDF_IGNORE, ".maxHaz",
				maxHazard_grid, tmp_grid, grid_cols, grid_rows, grid_cols_padded, Parptr->blx, Parptr->bly, Parptr->dx, output_params);
		}

	}
}
