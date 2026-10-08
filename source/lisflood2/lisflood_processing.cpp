/*
*****************************************************************************
ITERATEQ and UPDATEH
---------------------


*****************************************************************************
*/

#include "../lisflood.h"
#include "DataTypes.h"
#include "sgm_fast.h"
#include "file_tool.h"
#include "../utility.h"
#include <omp.h>

#include "../sgc.h"

#ifdef __unix__
#include <sched.h>
#include <numa.h>
#endif

void CopyToSubSubGridFlowInfo(int grid_index, int grid_index_next, int flow_index, NUMERIC_TYPE g,
	NUMERIC_TYPE grid_flow_distance, //grid_flow_distance is the distance between cell centres (meander coeffecient determines effective sg_flow_effective_distance)
	NUMERIC_TYPE grid_cell_width, // grid_cell_width cell width perpendicular to flow direction (For D8 dirns this is an 'effective width' to get the appropriate cell ratio)
	SubGridFlowInfo * flow_info, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr)
{
	int channel_group0 = Arrptr->SGCgroup[grid_index];
	int channel_group1 = Arrptr->SGCgroup[grid_index_next];
	NUMERIC_TYPE meander = C(0.5) * (SGCptr->SGCm[channel_group0] + SGCptr->SGCm[channel_group1]);
	flow_info->sg_flow_effective_distance[flow_index] = grid_flow_distance * meander;

	NUMERIC_TYPE g_friction_squared;
	if (Arrptr->SGCManningsn != NULL)
	{
		g_friction_squared = g * C(0.5) * Arrptr->SGCManningsn[grid_index] * Arrptr->SGCManningsn[grid_index_next];
	}
	else
	{
		g_friction_squared = g * SGCptr->SGCn[channel_group0];
	}

	flow_info->sg_flow_g_friction_sq[flow_index] = g_friction_squared;

	// PFU set constant channel ratio to minimum of two cells
	NUMERIC_TYPE width0 = Arrptr->SGCwidth[grid_index];
	NUMERIC_TYPE width1 = Arrptr->SGCwidth[grid_index_next];
	flow_info->sg_flow_ChannelRatio[flow_index] = getmin(width0,width1)/grid_cell_width;
	//printf("Debug channelratio: %d, %" NUM_FMT,flow_index,flow_info->sg_flow_ChannelRatio[flow_index]);
}

void CopyToSubGridCellInfo(int grid_cols_padded, int cell_index, int x, int y, SubGridCellInfo * sub_grid_cell_info, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr)
{
	int source_row_index = y * Parptr->xsz;
	int source_index = source_row_index + x;
	const NUMERIC_TYPE dx = Arrptr->dx[source_row_index];
	const NUMERIC_TYPE dy = Arrptr->dy[source_row_index];

	sub_grid_cell_info->sg_cell_x[cell_index] = x;
	sub_grid_cell_info->sg_cell_y[cell_index] = y;
	sub_grid_cell_info->sg_cell_grid_index_lookup[cell_index] = x + y * grid_cols_padded;
	
	sub_grid_cell_info->sg_cell_cell_area[cell_index] = Arrptr->dA[source_row_index];
	sub_grid_cell_info->sg_cell_dem[cell_index] = Arrptr->DEM[source_index];
	
	sub_grid_cell_info->sg_cell_SGC_width[cell_index] = Arrptr->SGCwidth[source_index];
	sub_grid_cell_info->sg_cell_SGC_BankFullHeight[cell_index] = Arrptr->SGCbfH[source_index];
	sub_grid_cell_info->sg_cell_SGC_BankFullVolume[cell_index] = Arrptr->SGCbfV[source_index];
	sub_grid_cell_info->sg_cell_SGC_c[cell_index] = Arrptr->SGCc[source_index];

	sub_grid_cell_info->sg_cell_SGC_group[cell_index] = Arrptr->SGCgroup[source_index];
	sub_grid_cell_info->sg_cell_SGC_is_large[cell_index] = (sub_grid_cell_info->sg_cell_SGC_width[cell_index] >= C(0.5)*(dx + dy));
}

int CheckSubGrid(const int source_index_this, const int i, const int j, const int grid_rows, const int grid_cols,
	const Arrays * Arrptr, const int weirs_enable, const int d8dirs_enable)
{
	int source_index_right = source_index_this + 1;
	int source_index_below = source_index_this + grid_cols;
	int source_index_belowright    = source_index_this + 1 + grid_cols;
	int source_index_belowleft     = source_index_this - 1 + grid_cols;

	int source_index_qx = (j * (grid_cols + 1)) + 1;
	int source_index_qy = ((j + 1) * (grid_cols + 1));

	int this_flow_count = 0;
	//printf("CheckSubGrid %d, %d\n", i, j);
	//fflush(NULL);
	// PFU use flow directions raster rather than connecting all adjacent cells
	// USES TauDEM directions: 1:E, 2:NE, 3:N, 4:NW, 5:W, 6:SW, 7:S, 8:SE
	// TODO, could code in option for multiple directions per cell. 
	// i.e. use arcgis dirns: 1:E, 2:SE, 4:S, 8:SW, 16:W, 32:NW, 64: N, 128: NE
	// Sum up values for multiple direction, 
	// Then check directions by checking the individual bits of the integer
	if (i < (grid_cols - 1) &&
		(Arrptr->SGCwidth[source_index_right] > C(0.0)) &&
		(Arrptr->DEM[source_index_right] != DEM_NO_DATA || Arrptr->ChanMask[source_index_right] > 0) && // typically the sub-grid model would not operate under DEM NoData, except if ChanMask is present
		(weirs_enable != ON || Arrptr->Weir_Identx[source_index_qx] == -1) ) // don't add the sub-grid calculation where there is a weir
	{
		//printf("sgx at %d, %d\n", i, j);
		// PFU, check direction raster (in both directions)
		if (Arrptr->SGCdirn == NULL)
		{
			this_flow_count++;
		}
		else
		{
			if (Arrptr->SGCdirn[source_index_this] == 1 || Arrptr->SGCdirn[source_index_right] == 5)
			{
				this_flow_count++;
			}
		}
	}
	if ((j < grid_rows - 1))
	{
		if ((Arrptr->SGCwidth[source_index_below] > C(0.0)) &&
			(Arrptr->DEM[source_index_below] != DEM_NO_DATA || Arrptr->ChanMask[source_index_below] > 0) &&
			(weirs_enable != ON || Arrptr->Weir_Identy[source_index_qy] == -1) ) // don't add the sub-grid calculation where there is a weir
		{
			//printf("sgy at %d, %d\n", i, j);
			if (Arrptr->SGCdirn == NULL)
			{
				this_flow_count++;
			}
			else
			{
				if(Arrptr->SGCdirn[source_index_this] == 7 || Arrptr->SGCdirn[source_index_below] == 3)
				{
					this_flow_count++;
				}
			}
		}
	}
	// PFU d8 dir belowright
	if (d8dirs_enable && (j < grid_rows - 1) && (i<grid_cols-1) )
	{
		if ((Arrptr->SGCwidth[source_index_belowright] > C(0.0)) &&
			(Arrptr->DEM[source_index_belowright] != DEM_NO_DATA || Arrptr->ChanMask[source_index_belowright] > 0)) 
		// warning: not checking for weirs in d8
		{
			//printf("sgy at %d, %d\n", i, j);
			if (Arrptr->SGCdirn == NULL)
			{
				this_flow_count++;
			}
			else
			{
				if (Arrptr->SGCdirn[source_index_this] == 8 || Arrptr->SGCdirn[source_index_belowright] == 4)
				{
					this_flow_count++;
				}
			}
		}
	}
	// PFU d8 dir belowleft
	if (d8dirs_enable && (j < grid_rows - 1) && (i>0) )
	{
		if ((Arrptr->SGCwidth[source_index_belowleft] > C(0.0)) &&
			(Arrptr->DEM[source_index_belowleft] != DEM_NO_DATA || Arrptr->ChanMask[source_index_belowleft] > 0) )
		// warning: not checking for weirs in d8
		{
			//printf("sgy at %d, %d\n", i, j);
			if (Arrptr->SGCdirn == NULL)
			{
				this_flow_count++;
			}
			else
			{
				if (Arrptr->SGCdirn[source_index_this] == 6 || Arrptr->SGCdirn[source_index_belowleft] == 2)
				{
					this_flow_count++;
				}
			}
		}
	}

	return this_flow_count;
}

///
/// sub grid divided into blocks
/// * allows rows to be processed in same open mp loop as flood plain
/// * not enough sub grid cells per row to fully utalize vectorization
///
void InitSubGridStructureByBlocks(SubGridState * sub_grid_state, SubGridRowList * sub_grid_layout,
	const int block_count,
	const int grid_cols_padded,
	const NUMERIC_TYPE g,
	int * flow_Qx_lookup_grid_tmp, int * flow_Qy_lookup_grid_tmp,
	NUMERIC_TYPE * Fp_xwidth,	NUMERIC_TYPE * Fp_ywidth,
	States *Statesptr, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr, const int verbose)
{
	if (verbose == ON)	printf("Init Sub Grid...\t");

	SubGridFlowInfo * sub_grid_flow_info = &sub_grid_layout->flow_info;
	SubGridCellInfo * sub_grid_cell_info = &sub_grid_layout->cell_info;

	const int grid_rows = Parptr->ysz;
	const int grid_cols = Parptr->xsz;

	//temporary data for building sub-grid structure
	int* sub_cell_lookup_grid_tmp = (int*)memory_allocate(sizeof(int)* grid_rows * Parptr->xsz);
	SetArrayValue(sub_cell_lookup_grid_tmp, -1, grid_rows * Parptr->xsz);

	sub_grid_layout->flow_row_count = (int*)memory_allocate(sizeof(int) * block_count);
	sub_grid_layout->cell_row_count = (int*)memory_allocate(sizeof(int) * block_count);

	int flow_count = 0;
	int cell_count = 0;

	//first count the channels 
	// the total amount of sub grid cells and sub grid flows must be know to divide them into blocks
	for (int j = 0; j < grid_rows; j++)
	{
		int flow_count_row = 0;
		int cell_count_row = 0;

		const int source_row_index = j * grid_cols;
		// Intel compiler 2015 causes a segfault with this loop with simd vector
#pragma novector
		//#pragma ivdep
		//#pragma simd reduction(+:flow_count_row, cell_count_row)
		for (int i = 0; i < grid_cols; i++)
		{
			int source_index_this = source_row_index + i;

			if (Arrptr->SGCwidth[source_index_this] > C(0.0) && (Arrptr->DEM[source_index_this] != DEM_NO_DATA || Arrptr->ChanMask[source_index_this] > 0))
			{
				int this_flow_count;
				this_flow_count = CheckSubGrid(source_index_this, i, j, grid_rows, grid_cols, Arrptr, Statesptr->weirs, Statesptr->SGCd8);
				flow_count_row += this_flow_count;
				cell_count_row++;
			}
		}
		cell_count += cell_count_row;
		flow_count += flow_count_row;
	}

	float target_flows_per_block = (float)flow_count / block_count;
	//int cells_per_block = (int)ceil((double)cell_count / block_count);

	/// each block will take flows from a range of rows
	IndexRange * row_range = (IndexRange*)memory_allocate(block_count * sizeof(IndexRange));
	int* block_count_flows = (int*)memory_allocate(block_count * sizeof(int));
	int* block_count_cells = (int*)memory_allocate(block_count * sizeof(int));
	memset(block_count_flows, 0, block_count * sizeof(int));
	memset(block_count_cells, 0, block_count * sizeof(int));
	for (int i = 0; i < block_count; i++)
	{
		row_range[i].start = -1;
		row_range[i].end = -1;
	}
	int block_index = 0;
	int allocated_flows = 0;

	if (flow_count != 0)
	{
		row_range[0].start = 0;
		for (int j = 0; j < grid_rows; j++)
		{
			// increment blocks only on row boundary
			if (block_index < (block_count - 1) &&
				allocated_flows >(block_index + 1) * target_flows_per_block)
			{
				row_range[block_index].end = j + 1;
				block_index++;
				row_range[block_index].start = j + 1;
			}

			const int source_row_index = j * grid_cols;
			// Intel compiler 2015 causes a segfault with this loop with simd vector
#pragma novector
			//#pragma ivdep
			//#pragma simd reduction(+:flow_count_row, cell_count_row)
			for (int i = 0; i < grid_cols; i++)
			{
				int source_index_this = source_row_index + i;

				if (Arrptr->SGCwidth[source_index_this] > C(0.0) && (Arrptr->DEM[source_index_this] != DEM_NO_DATA || Arrptr->ChanMask[source_index_this] > 0))
				{
					int this_flow_count;
					this_flow_count = CheckSubGrid(source_index_this, i, j, grid_rows, grid_cols, Arrptr, Statesptr->weirs, Statesptr->SGCd8);
					allocated_flows += this_flow_count;
					block_count_flows[block_index] += this_flow_count;
					block_count_cells[block_index]++;
				}
			}
		}

		row_range[block_index].end = grid_rows;
	}


#ifdef _DEBUG
	printf("\n");
	for (int bi = 0; bi < block_count; bi++)
	{
		int rows_in_block = row_range[bi].end - row_range[bi].start;
		NUMERIC_TYPE percent_rows_in_block = ((NUMERIC_TYPE)rows_in_block / grid_rows) * 100;

		int flows_this_block = block_count_flows[bi];
		int cells_this_block = block_count_cells[bi];

		NUMERIC_TYPE percent_work1 = flow_count > 0 ? ((NUMERIC_TYPE)flows_this_block / flow_count) * 100 : 0;
		NUMERIC_TYPE percent_work2 = cell_count > 0 ? ((NUMERIC_TYPE)cells_this_block / cell_count) * 100 : 0;
		printf("block: %d row: %d -> %d [%d](%.2" NUM_FMT"%%) work: [%d : %.2" NUM_FMT"%%] [%d : %.2" NUM_FMT"%%]\n",
			bi,
			row_range[bi].start, row_range[bi].end,
			rows_in_block, percent_rows_in_block, flows_this_block, percent_work1, cells_this_block, percent_work2);

	}
#endif

	int block_max_flows = 0;
	int block_max_cells = 0;
	for (int bi = 0; bi < block_count; bi++)
	{
		block_max_flows = max(block_count_flows[bi], block_max_flows);
		block_max_cells = max(block_count_cells[bi], block_max_cells);
	}

	int block_memory_size = max(block_max_flows, block_max_cells);
	block_memory_size = block_memory_size + (GRID_ALIGN_WIDTH - (block_memory_size % GRID_ALIGN_WIDTH)) % GRID_ALIGN_WIDTH;
	int memory_size = block_memory_size * block_count;
	int flow_memory_size = memory_size;
	int cell_memory_size = memory_size;

	//sub_grid_flow_info->flow_count = flow_count;
	sub_grid_layout->row_cols_padded = block_memory_size;

	sub_grid_state->sg_flow_Q = (NUMERIC_TYPE*)memory_allocate(flow_memory_size * sizeof(NUMERIC_TYPE));
	sub_grid_flow_info->sg_flow_ChannelRatio = (NUMERIC_TYPE*)memory_allocate(flow_memory_size * sizeof(NUMERIC_TYPE));
	sub_grid_flow_info->sg_flow_effective_distance = (NUMERIC_TYPE*)memory_allocate(flow_memory_size * sizeof(NUMERIC_TYPE));
	sub_grid_flow_info->sg_flow_g_friction_sq = (NUMERIC_TYPE*)memory_allocate(flow_memory_size * sizeof(NUMERIC_TYPE));
	sub_grid_flow_info->sg_pair_cell_index_lookup = (int*)memory_allocate((2 * block_memory_size) * block_count * sizeof(int));

	AllocateSubGridCellInfo((2 * block_memory_size) * grid_rows, &sub_grid_flow_info->flow_pair);

	sub_grid_flow_info->sg_cell_flow_lookup = (SubGridFlowLookup*)memory_allocate(cell_memory_size * sizeof(SubGridFlowLookup));
	AllocateSubGridCellInfo(cell_memory_size, sub_grid_cell_info);
	sub_grid_cell_info->cell_count = cell_count;

	//initialise in thread - to set the affinity via first touch
#pragma omp parallel for default(shared) schedule(static)
	for (int bi = 0; bi < block_count; bi++)
	{
		int sg_block_start_index = bi * block_memory_size;
		int flow_count_block = 0;
		int cell_count_block = 0;

		int cell_index = sg_block_start_index;
		for (int j = row_range[bi].start; j < row_range[bi].end; j++)
		{
			const int source_row_index = j * Parptr->xsz;
			const int padded_grid_row_index = j * grid_cols_padded;
			const bool not_last_row = (j < Parptr->ysz - 1);


			for (int i = 0; i < Parptr->xsz; i++)
			{
				int source_index_this = source_row_index + i;

				if (Arrptr->SGCwidth[source_index_this] > C(0.0) && (Arrptr->DEM[source_index_this] != DEM_NO_DATA || Arrptr->ChanMask[source_index_this] > 0))
				{
					int this_flow_count = CheckSubGrid(source_index_this, i, j, grid_rows, Parptr->xsz, Arrptr, Statesptr->weirs, Statesptr->SGCd8);
					flow_count_block += this_flow_count;
					cell_count_block++;
					sub_cell_lookup_grid_tmp[source_index_this] = cell_index;

					//these will be used for d8 diagonal flows
					sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_add[0] = -1;
					sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_add[1] = -1;
					sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_subtract[0] = -1;
					sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_subtract[1] = -1;

					CopyToSubGridCellInfo(grid_cols_padded, cell_index, i, j, sub_grid_cell_info, Parptr, Arrptr, SGCptr);

					cell_index++;
				}
				else
				{
					//Arrptr->SGCbfH[source_index] = 0;
				}
			}
		}
		sub_grid_layout->flow_row_count[bi] = flow_count_block;
		sub_grid_layout->cell_row_count[bi] = cell_count_block;
	}
	// second loop to ensure sub_cell_lookup_grid_tmp fully initialized before using it.
	// initialise in thread - to set the affinity via first touch
#pragma omp parallel for default(shared) schedule(static)
	for (int bi = 0; bi < block_count; bi++)
	{
		int sg_block_start_index = bi * block_memory_size;
		int sg_block_pair_start_index = bi * 2 * block_memory_size;
		int row_flow_index = 0;

		int flow_count_block = 0;
		int cell_count_block = 0;

		for (int j = row_range[bi].start; j < row_range[bi].end; j++)
		{
			const int source_row_index = j * Parptr->xsz;
			const bool not_last_row = (j < Parptr->ysz - 1);
			const NUMERIC_TYPE row_dx = Arrptr->dx[source_row_index];
			const NUMERIC_TYPE row_dy = Arrptr->dy[source_row_index];

			for (int i = 0; i < Parptr->xsz; i++)
			{
				int source_index_this = source_row_index + i;
				int source_index_right = source_index_this + 1;
				int source_index_below = source_index_this + Parptr->xsz;

				int source_index_qx = (j * (grid_cols + 1)) + 1;
				int source_index_qy = ((j + 1) * (grid_cols + 1));

				if (Arrptr->SGCwidth[source_index_this] > C(0.0) && (Arrptr->DEM[source_index_this] != DEM_NO_DATA || Arrptr->ChanMask[source_index_this] > 0))
				{
					if (i < (grid_cols - 1) &&
						(Arrptr->SGCwidth[source_index_right] > C(0.0)) &&
						(Arrptr->DEM[source_index_right] != DEM_NO_DATA || Arrptr->ChanMask[source_index_right] > 0) &&
						(Statesptr->weirs != ON || Arrptr->Weir_Identx[source_index_qx] == -1)) // don't add the sub-grid calculation where there is a weir)
					{
						int flow_index = sg_block_start_index + row_flow_index;
						int flow_pair_index = sg_block_pair_start_index + 2 * row_flow_index;

						int grid_index = j * grid_cols_padded + i;
						//store the index of the q by grid, weirs/bridges need to know the upstream flow
						flow_Qx_lookup_grid_tmp[grid_index] = flow_index;

						sub_grid_state->sg_flow_Q[flow_index] = C(0.0);
						sub_grid_flow_info->sg_flow_effective_distance[flow_index] = C(0.0);
						sub_grid_flow_info->sg_flow_g_friction_sq[flow_index] = C(0.0);

						ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index);
						ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index + 1);

						// Qx
						CopyToSubSubGridFlowInfo(source_index_this, source_index_right, flow_index, g, row_dx, row_dy, sub_grid_flow_info, Parptr, Arrptr, SGCptr);
						CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index, i, j, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);
						CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index + 1, i + 1, j, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);

						int cell_list_index_this = sub_cell_lookup_grid_tmp[source_index_this];
						int cell_list_index_right = sub_cell_lookup_grid_tmp[source_index_right];

						sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index] = cell_list_index_this;
						sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index + 1] = cell_list_index_right;

						sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_this].flow_subtract[0] = flow_index;//subtract from this cell
						sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_right].flow_add[0] = flow_index;//add to the right

						row_flow_index++;
					}
					if (not_last_row)
					{
						if ((Arrptr->SGCwidth[source_index_below] > C(0.0)) &&
							(Arrptr->DEM[source_index_below] != DEM_NO_DATA || Arrptr->ChanMask[source_index_below] > 0) &&
							(Statesptr->weirs != ON || Arrptr->Weir_Identy[source_index_qy] == -1)) // don't add the sub-grid calculation where there is a weir)
						{
							int flow_index = sg_block_start_index + row_flow_index;
							int flow_pair_index = sg_block_pair_start_index + 2 * row_flow_index;

							//store the index of the q by grid, weirs/bridges need to know the upstream flow
							int grid_index = j * grid_cols_padded + i;
							flow_Qy_lookup_grid_tmp[grid_index] = flow_index;

							sub_grid_state->sg_flow_Q[flow_index] = C(0.0);
							sub_grid_flow_info->sg_flow_effective_distance[flow_index] = C(0.0);
							sub_grid_flow_info->sg_flow_g_friction_sq[flow_index] = C(0.0);

							ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index);
							ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index + 1);

							// Qy
							CopyToSubSubGridFlowInfo(source_index_this, source_index_below, flow_index, g, row_dy, row_dx, sub_grid_flow_info, Parptr, Arrptr, SGCptr);
							CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index, i, j, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);
							CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index + 1, i, j + 1, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);

							int cell_list_index_this = sub_cell_lookup_grid_tmp[source_index_this];
							int cell_list_index_below = sub_cell_lookup_grid_tmp[source_index_below];

							sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index] = cell_list_index_this;
							sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index + 1] = cell_list_index_below;

							sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_this].flow_subtract[1] = flow_index; //subtract from this cell
							sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_below].flow_add[1] = flow_index; //add to the cell below

							row_flow_index++;
						}
					}
				}
			}
		}
	}




	if (verbose == ON)	printf("sub grid flows: %d sub grid cells:%d\n", flow_count, cell_count);

#ifdef _DEBUG

	for (int bi = 0; bi < block_count; bi++)
	{
		printf("block %d flows: %d cells %d\n", bi, sub_grid_layout->flow_row_count[bi], sub_grid_layout->cell_row_count[bi]);
	}

#endif

	memory_free(&sub_cell_lookup_grid_tmp);
}

///
/// sub grid divided into rows
/// * allows rows to be processed in same open mp loop as flood plain
/// * not enough sub grid cells per row to fully utalize vectorization
///
void InitSubGridStructureByRows(SubGridState * sub_grid_state, SubGridRowList * sub_grid_layout,
	const int grid_cols_padded,
	const NUMERIC_TYPE g,
	int * flow_Qx_lookup_grid_tmp, int * flow_Qy_lookup_grid_tmp,
	NUMERIC_TYPE * Fp_xwidth,	NUMERIC_TYPE * Fp_ywidth,
	States *Statesptr, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr, const int verbose)
{
	if (verbose == ON)	printf("Init Sub Grid...\t");

	SubGridFlowInfo * sub_grid_flow_info = &sub_grid_layout->flow_info;
	SubGridCellInfo * sub_grid_cell_info = &sub_grid_layout->cell_info;

	const int grid_rows = Parptr->ysz;
	const int grid_cols = Parptr->xsz;

	//temporary data for building sub-grid structure
	int* sub_cell_lookup_grid_tmp = (int*)memory_allocate(sizeof(int)* grid_rows * Parptr->xsz);
	SetArrayValue(sub_cell_lookup_grid_tmp, -1, grid_rows * Parptr->xsz);

	sub_grid_layout->flow_row_count = (int*)memory_allocate(sizeof(int) * grid_rows);
	sub_grid_layout->cell_row_count = (int*)memory_allocate(sizeof(int) * grid_rows);

	int flow_count = 0;
	int cell_count = 0;
	int flow_count_row_max = 0;
	int cell_count_row_max = 0;


	//first count the channels - to know how much memory to allocate
	for (int j = 0; j < grid_rows; j++)
	{
		int flow_count_row = 0;
		int cell_count_row = 0;

		const int source_row_index = j * grid_cols;


		//PFU initialise fp widths to 1
		// Y-grid has one extra row (for xwidths)
		for (int i = 0; i < grid_cols_padded; i++){
			if (j==0)
				Fp_xwidth[i] = C(1.0);
			Fp_ywidth[i + j * grid_cols_padded] = C(1.0);
			Fp_xwidth[i + (j+1) * grid_cols_padded] = C(1.0);
		}


		// Intel compiler 2015 causes a segfault with this loop with simd vector
#pragma novector
		//#pragma ivdep
		//#pragma simd reduction(+:flow_count_row, cell_count_row)
		for (int i = 0; i < grid_cols; i++)
		{
			int source_index_this = source_row_index + i;

			if (Arrptr->SGCwidth[source_index_this] > C(0.0) && (Arrptr->DEM[source_index_this] != DEM_NO_DATA || Arrptr->ChanMask[source_index_this] > 0))
			{
				int this_flow_count;
				this_flow_count = CheckSubGrid(source_index_this, i, j, grid_rows, grid_cols, Arrptr, Statesptr->weirs, Statesptr->SGCd8);
				flow_count_row += this_flow_count;
				cell_count_row++;
			}

		}
		sub_grid_layout->flow_row_count[j] = flow_count_row;
		sub_grid_layout->cell_row_count[j] = cell_count_row;

		if (cell_count_row > cell_count_row_max)
			cell_count_row_max = cell_count_row;
		if (flow_count_row > flow_count_row_max)
			flow_count_row_max = flow_count_row;

		cell_count += cell_count_row;
		flow_count += flow_count_row;
	}

	int row_cols_padded = max(flow_count_row_max, cell_count_row_max);
	row_cols_padded = row_cols_padded + (GRID_ALIGN_WIDTH - (row_cols_padded % GRID_ALIGN_WIDTH)) % GRID_ALIGN_WIDTH;
	int memory_size = row_cols_padded * grid_rows;
	int flow_memory_size = memory_size;
	int cell_memory_size = memory_size;

	//sub_grid_flow_info->flow_count = flow_count;
	sub_grid_layout->row_cols_padded = row_cols_padded;

	sub_grid_state->sg_flow_Q = (NUMERIC_TYPE*)memory_allocate(flow_memory_size * sizeof(NUMERIC_TYPE));
	sub_grid_flow_info->sg_flow_ChannelRatio = (NUMERIC_TYPE*)memory_allocate(flow_memory_size * sizeof(NUMERIC_TYPE));
	sub_grid_state->sg_velocity = (NUMERIC_TYPE*)memory_allocate(flow_memory_size * sizeof(NUMERIC_TYPE));

	sub_grid_flow_info->sg_flow_effective_distance = (NUMERIC_TYPE*)memory_allocate(flow_memory_size * sizeof(NUMERIC_TYPE));
	sub_grid_flow_info->sg_flow_g_friction_sq = (NUMERIC_TYPE*)memory_allocate(flow_memory_size * sizeof(NUMERIC_TYPE));
	sub_grid_flow_info->sg_pair_cell_index_lookup = (int*)memory_allocate((2 * row_cols_padded) * grid_rows * sizeof(int));

	AllocateSubGridCellInfo((2 * row_cols_padded) * grid_rows, &sub_grid_flow_info->flow_pair);

	sub_grid_flow_info->sg_cell_flow_lookup = (SubGridFlowLookup*)memory_allocate(cell_memory_size * sizeof(SubGridFlowLookup));
	AllocateSubGridCellInfo(cell_memory_size, sub_grid_cell_info);
	sub_grid_cell_info->cell_count = cell_count;

	//initialise in thread - to set the affinity via first touch
#pragma omp parallel for default(shared) schedule(static)
	for (int j = 0; j < grid_rows; j++)
	{
		int sg_row_start_index = j * row_cols_padded;

		const int source_row_index = j * Parptr->xsz;
		const int padded_grid_row_index = j * grid_cols_padded;
		const bool not_last_row = (j < Parptr->ysz - 1);

		int cell_count_row = 0;
		int flow_count_row = 0;

		int cell_index = sg_row_start_index;
		for (int i = 0; i < Parptr->xsz; i++)
		{
			int source_index_this = source_row_index + i;

			if (Arrptr->SGCwidth[source_index_this] > C(0.0) && (Arrptr->DEM[source_index_this] != DEM_NO_DATA || Arrptr->ChanMask[source_index_this] > 0))
			{
				int this_flow_count = CheckSubGrid(source_index_this, i, j, grid_rows, Parptr->xsz, Arrptr, Statesptr->weirs, Statesptr->SGCd8);
				flow_count_row += this_flow_count;
				cell_count_row++;
				sub_cell_lookup_grid_tmp[source_index_this] = cell_index;

				//these will be used for d8 diagonal flows
				sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_add[0] = -1;
				sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_add[1] = -1;
				sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_subtract[0] = -1;
				sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_subtract[1] = -1;
				// PFU initialize D8 flow add and subtract
				// For now do this also for D4 case (otherwise we need extra checks in sgm_fast for SGCd8 flag)
				//if (Statesptr->SGCd8) // PFU add d8 directions (BUT are they used at all?)
				//{
				sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_add[2] = -1;
				sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_add[3] = -1;
				sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_subtract[2] = -1;
				sub_grid_flow_info->sg_cell_flow_lookup[cell_index].flow_subtract[3] = -1;
				//}

				CopyToSubGridCellInfo(grid_cols_padded, cell_index, i, j, sub_grid_cell_info, Parptr, Arrptr, SGCptr);
#ifdef _DEBUG
				if (sub_grid_cell_info->sg_cell_SGC_is_large[cell_index])
					printf("large sub grid %d %d\n", i, j);
#endif
				cell_index++;
			}
			else
			{
				//Arrptr->SGCbfH[source_index] = 0;
			}
		}
	}
	// second loop to ensure sub_cell_lookup_grid_tmp fully initialized before using it.
	// initialise in thread - to set the affinity via first touch
#pragma omp parallel for default(shared) schedule(static)
	for (int j = 0; j < grid_rows; j++)
	{
		const int source_row_index = j * Parptr->xsz;
		const bool not_last_row = (j < Parptr->ysz - 1);
		const NUMERIC_TYPE row_dx = Arrptr->dx[source_row_index];
		const NUMERIC_TYPE row_dy = Arrptr->dy[source_row_index];

		int sg_row_start_index = j * row_cols_padded;
		int sg_row_pair_start_index = j * 2 * row_cols_padded;
		int row_flow_index = 0;

		for (int i = 0; i < Parptr->xsz; i++)
		{
			int source_index_this = source_row_index + i;
			int source_index_right = source_index_this + 1;
			int source_index_below = source_index_this + Parptr->xsz;
			int source_index_belowright = source_index_this + 1 + Parptr->xsz;
			int source_index_belowleft = source_index_this - 1 + Parptr->xsz;

			int source_index_right_q = (j * (Parptr->xsz + 1)) + 1;
			int source_index_below_q = ((j + 1) * (Parptr->xsz + 1));

			if (Arrptr->SGCwidth[source_index_this] > C(0.0) && (Arrptr->DEM[source_index_this] != DEM_NO_DATA || Arrptr->ChanMask[source_index_this] > 0))
			{
				// Link with cell to the right
				if (i < (Parptr->xsz - 1) &&
					(Arrptr->SGCwidth[source_index_right] > C(0.0)) &&
					(Arrptr->DEM[source_index_right] != DEM_NO_DATA || Arrptr->ChanMask[source_index_right] > 0)  &&
					(Statesptr->weirs != ON || Arrptr->Weir_Identx[source_index_right_q] == -1)) // don't add the sub-grid calculation where there is a weir)
				{
					// check directions raster for link
					bool dirn_allowed = 1;
					if (Arrptr->SGCdirn != NULL)
					{
						if  (Arrptr->SGCdirn[source_index_this] != 1 && Arrptr->SGCdirn[source_index_right] != 5 ) dirn_allowed = 0;
					}
					if (dirn_allowed)
					{
						int flow_index = sg_row_start_index + row_flow_index;
						int flow_pair_index = sg_row_pair_start_index + 2 * row_flow_index;
						// printf("Link right: %i-%i\n",source_index_this,source_index_right);

						int grid_index = j * grid_cols_padded + i;
						//store the index of the q by grid, weirs/bridges need to know the upstream flow
						flow_Qx_lookup_grid_tmp[grid_index] = flow_index;

						sub_grid_state->sg_flow_Q[flow_index] = C(0.0);
						sub_grid_flow_info->sg_flow_effective_distance[flow_index] = C(0.0);
						sub_grid_flow_info->sg_flow_g_friction_sq[flow_index] = C(0.0);

						ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index);
						ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index + 1);

						// Qx
						CopyToSubSubGridFlowInfo(source_index_this, source_index_right, flow_index, g, row_dx, row_dy, sub_grid_flow_info, Parptr, Arrptr, SGCptr);
						CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index, i, j, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);
						CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index + 1, i + 1, j, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);

						int cell_list_index_this = sub_cell_lookup_grid_tmp[source_index_this];
						int cell_list_index_right = sub_cell_lookup_grid_tmp[source_index_right];

						sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index] = cell_list_index_this;
						sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index + 1] = cell_list_index_right;

						sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_this].flow_subtract[0] = flow_index;//subtract from this cell
						sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_right].flow_add[0] = flow_index;//add to the right

						// PFU subtract channel ratio from floodplain width (floodplain width initialized to 1)
						// Widths are on Q grid (between floodplain cells)
						// use omp critical to prevent different threads writing to same width variable
						#pragma omp critical(widthy)
						{
						Fp_ywidth[grid_index+1] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index];
						}

						row_flow_index++;
					}
				}
				// Link with cell below
				if (not_last_row)
				{
					if ((Arrptr->SGCwidth[source_index_below] > C(0.0)) &&
						(Arrptr->DEM[source_index_below] != DEM_NO_DATA || Arrptr->ChanMask[source_index_below] > 0) &&
						(Statesptr->weirs != ON || Arrptr->Weir_Identy[source_index_below_q] == -1)) // don't add the sub-grid calculation where there is a weir)
					{
						// check directions raster for link
						bool dirn_allowed = 1;
						if (Arrptr->SGCdirn != NULL)
						{
							if ( Arrptr->SGCdirn[source_index_this] != 7 && Arrptr->SGCdirn[source_index_right] != 3 ) dirn_allowed = 0;
						}
						if (dirn_allowed)
							{

							int flow_index = sg_row_start_index + row_flow_index;
							int flow_pair_index = sg_row_pair_start_index + 2 * row_flow_index;
							// printf("Link below: %i-%i\n",source_index_this,source_index_below);

							//store the index of the q by grid, weirs/bridges need to know the upstream flow
							int grid_index = j * grid_cols_padded + i;
							flow_Qy_lookup_grid_tmp[grid_index] = flow_index;

							sub_grid_state->sg_flow_Q[flow_index] = C(0.0);
							sub_grid_flow_info->sg_flow_effective_distance[flow_index] = C(0.0);
							sub_grid_flow_info->sg_flow_g_friction_sq[flow_index] = C(0.0);

							ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index);
							ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index + 1);

							// Qy
							CopyToSubSubGridFlowInfo(source_index_this, source_index_below, flow_index, g, row_dy, row_dx, sub_grid_flow_info, Parptr, Arrptr, SGCptr);
							CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index, i, j, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);
							CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index + 1, i, j + 1, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);

							int cell_list_index_this = sub_cell_lookup_grid_tmp[source_index_this];
							int cell_list_index_below = sub_cell_lookup_grid_tmp[source_index_below];

							sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index] = cell_list_index_this;
							sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index + 1] = cell_list_index_below;

							sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_this].flow_subtract[1] = flow_index; //subtract from this cell
							sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_below].flow_add[1] = flow_index; //add to the cell below

							// PFU subtract channel ratio from floodplain width (floodplain width initialized to 1)
							// Widths are on Q grid (between floodplain cells)
							// use omp critical to prevent different threads writing to same width variable
							#pragma omp critical(widthx)
							{
							Fp_xwidth[grid_index+grid_cols_padded] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index];
							}

							row_flow_index++;
						}
					}
				}
				// Link with the cell belowright (PFU d8 directions)
				if (Statesptr->SGCd8 && not_last_row && (i < (Parptr->xsz - 1)))
				{
					if ((Arrptr->SGCwidth[source_index_belowright] > C(0.0)) &&
						(Arrptr->DEM[source_index_belowright] != DEM_NO_DATA || Arrptr->ChanMask[source_index_belowright] > 0)) // PFU: check directions raster for link
					{
						// check directions raster for link
						bool dirn_allowed = 1;
						if (Arrptr->SGCdirn != NULL)
						{
							if ( Arrptr->SGCdirn[source_index_this] != 8 && Arrptr->SGCdirn[source_index_right] != 4 ) dirn_allowed = 0;
						}
						if (dirn_allowed)
						{

							int flow_index = sg_row_start_index + row_flow_index;
							int flow_pair_index = sg_row_pair_start_index + 2 * row_flow_index;
							// printf("Link belowright: %i-%i\n",source_index_this,source_index_belowright);
							//store the index of the q by grid, weirs/bridges need to know the upstream flow
							int grid_index = j * grid_cols_padded + i;
							flow_Qy_lookup_grid_tmp[grid_index] = flow_index;

							sub_grid_state->sg_flow_Q[flow_index] = C(0.0);
							sub_grid_flow_info->sg_flow_ChannelRatio[flow_index] = C(0.0);
							sub_grid_flow_info->sg_flow_effective_distance[flow_index] = C(0.0);
							sub_grid_flow_info->sg_flow_g_friction_sq[flow_index] = C(0.0);

							ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index);
							ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index + 1);

							NUMERIC_TYPE flow_dist = SQRT(row_dy*row_dy + row_dx*row_dx);
							// This cell width (perpendicular to the flow) is used for the channel_ratio 
							// Redefine 'effective' cell width for diagonal flows, to get appropriate channel ratio (channel ratio is channel_width / cell_width)). 
							NUMERIC_TYPE cell_width = (row_dy * row_dx) / flow_dist; 
							CopyToSubSubGridFlowInfo(source_index_this, source_index_belowright, flow_index, g, flow_dist, cell_width, sub_grid_flow_info, Parptr, Arrptr, SGCptr);
							CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index, i, j, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);
							CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index + 1, i + 1, j + 1, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);

							int cell_list_index_this = sub_cell_lookup_grid_tmp[source_index_this];
							int cell_list_index_belowright = sub_cell_lookup_grid_tmp[source_index_belowright];

							sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index] = cell_list_index_this;
							sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index + 1] = cell_list_index_belowright;

							sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_this].flow_subtract[2] = flow_index; //subtract from this cell
							sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_belowright].flow_add[2] = flow_index; //add to the cell belowright

							// This is the ratio of dx/dy cell sides to the diagonal cell length
							NUMERIC_TYPE dx_diag_ratio = row_dx/flow_dist;
							NUMERIC_TYPE dy_diag_ratio = row_dy/flow_dist;

							// PFU subtract channel ratio from floodplain width (floodplain width initialized to 1)
							// Widths are on Q grid (between floodplain cells)
							// use omp critical to prevent different threads writing to same width variable
							#pragma omp critical(widthx)
							{
							Fp_xwidth[grid_index+grid_cols_padded] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index]/C(2.0)*dx_diag_ratio;
							Fp_xwidth[grid_index+grid_cols_padded+1] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index]/C(2.0)*dx_diag_ratio;
							}
							#pragma omp critical(widthy)
							{
							Fp_ywidth[grid_index+1] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index]/C(2.0)*dy_diag_ratio;
							Fp_ywidth[grid_index+1+grid_cols_padded] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index]/C(2.0)*dy_diag_ratio;
							}

							row_flow_index++;
						}
					}
				}
				// Link with the cell belowleft (PFU d8 directions)
				if (Statesptr->SGCd8 && not_last_row && (i >0))
				{
					if ((Arrptr->SGCwidth[source_index_belowleft] > C(0.0)) &&
						(Arrptr->DEM[source_index_belowleft] != DEM_NO_DATA || Arrptr->ChanMask[source_index_belowleft] > 0))
					{
						// check directions raster for link
						bool dirn_allowed = 1;
						if (Arrptr->SGCdirn != NULL)
						{
							if ( Arrptr->SGCdirn[source_index_this] != 6 && Arrptr->SGCdirn[source_index_right] != 2 ) dirn_allowed = 0;
						}
						if (dirn_allowed)
						{

							int flow_index = sg_row_start_index + row_flow_index;
							int flow_pair_index = sg_row_pair_start_index + 2 * row_flow_index;
							// printf("Link belowleft: %i-%i\n",source_index_this,source_index_belowleft);

							//store the index of the q by grid, weirs/bridges need to know the upstream flow
							int grid_index = j * grid_cols_padded + i;
							flow_Qy_lookup_grid_tmp[grid_index] = flow_index;

							sub_grid_state->sg_flow_Q[flow_index] = C(0.0);
							sub_grid_flow_info->sg_flow_effective_distance[flow_index] = C(0.0);
							sub_grid_flow_info->sg_flow_g_friction_sq[flow_index] = C(0.0);

							ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index);
							ZeroSubGridCellInfo(&sub_grid_flow_info->flow_pair, flow_pair_index + 1);

							NUMERIC_TYPE flow_dist = SQRT(row_dy*row_dy + row_dx*row_dx);
							// This cell width (perpendicular to the flow) is used for the channel_ratio 
							// Redefine 'effective' cell width for diagonal flows, to get appropriate channel ratio (channel ratio is channel_width / cell_width)). 
							NUMERIC_TYPE cell_width = (row_dy * row_dx) / flow_dist; 
							CopyToSubSubGridFlowInfo(source_index_this, source_index_belowleft, flow_index, g, flow_dist, cell_width, sub_grid_flow_info, Parptr, Arrptr, SGCptr);
							CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index, i, j, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);
							CopyToSubGridCellInfo(grid_cols_padded, flow_pair_index + 1, i - 1, j + 1, &sub_grid_flow_info->flow_pair, Parptr, Arrptr, SGCptr);

							int cell_list_index_this = sub_cell_lookup_grid_tmp[source_index_this];
							int cell_list_index_belowleft = sub_cell_lookup_grid_tmp[source_index_belowleft];

							sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index] = cell_list_index_this;
							sub_grid_flow_info->sg_pair_cell_index_lookup[flow_pair_index + 1] = cell_list_index_belowleft;

							sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_this].flow_subtract[3] = flow_index; //subtract from this cell
							sub_grid_flow_info->sg_cell_flow_lookup[cell_list_index_belowleft].flow_add[3] = flow_index; //add to the cell belowleft


							// This is the ratio of dx/dy cell sides to the diagonal cell length
							NUMERIC_TYPE dx_diag_ratio = row_dx/flow_dist;
							NUMERIC_TYPE dy_diag_ratio = row_dy/flow_dist;

							// PFU subtract channel ratio from floodplain width (floodplain width initialized to 1)
							// Widths are on Q grid (between floodplain cells)
							// use omp critical to prevent different threads writing to same width variable
							#pragma omp critical(widthx)
							{
							Fp_xwidth[grid_index+grid_cols_padded-1] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index]/C(2.0)*dx_diag_ratio;
							Fp_xwidth[grid_index+grid_cols_padded] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index]/C(2.0)*dx_diag_ratio;
							}
							#pragma omp critical(widthy)
							{
							Fp_ywidth[grid_index] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index]/C(2.0)*dy_diag_ratio;
							Fp_ywidth[grid_index+grid_cols_padded] -= sub_grid_flow_info->sg_flow_ChannelRatio[flow_index]/C(2.0)*dy_diag_ratio;
							}

							row_flow_index++;
						}
					}
				}
			}
		}
	}

	#pragma omp parallel for default(shared) schedule(static)
	for (int j = 0; j < grid_rows; j++)
	{
		const int padded_grid_row_index = j * grid_cols_padded;
		const int padded_grid_nextrow_index = (j+1) * grid_cols_padded;	
		const int source_row_index = j * Parptr->xsz;
		const NUMERIC_TYPE row_dx = Arrptr->dx[source_row_index];
		const NUMERIC_TYPE row_dy = Arrptr->dy[source_row_index];
		//PFU reset negative widths to 0, multiply by cell width
		// Y-grid has one extra row
		for (int i = 0; i < grid_cols_padded; i++){
			if (j==0)
				Fp_xwidth[i] = getmax(Fp_xwidth[i],C(0.0))*row_dx;
			//if (Fp_xwidth[i + padded_grid_row_index] < C(1.0))
			//	printf("xwidth: %" NUM_FMT,getmax(Fp_xwidth[i + padded_grid_row_index],C(0.0))*row_dx);
			Fp_ywidth[i + padded_grid_row_index] = getmax(Fp_ywidth[i + padded_grid_row_index],C(0.0))*row_dy;
			//if (Fp_ywidth[i + padded_grid_nextrow_index] < C(1.0))
			//	printf("ywidth: %" NUM_FMT,getmax(Fp_ywidth[i + padded_grid_nextrow_index],C(0.0))*row_dx);
			Fp_xwidth[i + padded_grid_nextrow_index] = getmax(Fp_xwidth[i + padded_grid_nextrow_index],C(0.0))*row_dx;

		}
	}

	if (verbose == ON)	printf("sub grid flows: %d sub grid cells:%d\n", flow_count, cell_count);

	memory_free(&sub_cell_lookup_grid_tmp);
}

void PointSourceToWaterSource(int ps_index, int ws_index, const int grid_cols_padded, const int grid_cols,
	const NUMERIC_TYPE g,
	WaterSource *ps_info,
	Pars *Parptr, Arrays *Arrptr,
	BoundCs * BCptr, SGCprams *SGCptr)
{
	ps_info->Ident[ws_index] = BCptr->PS_Ident[ps_index];
	ps_info->Val[ws_index] = BCptr->PS_Val[ps_index];
	ps_info->timeSeries[ws_index] = BCptr->PS_TimeSeries[ps_index];

	ps_info->Q_FP_old[ws_index] = BCptr->PS_Q_FP_old[ps_index];
	ps_info->Q_SG_old[ws_index] = BCptr->PS_Q_SG_old[ps_index];;

	int x = BCptr->xpi[ps_index];
	int y = BCptr->ypi[ps_index];

	// reading from non padded grids, so use grid_cols
	int source_index = y * grid_cols + x;

	NUMERIC_TYPE n;
	if (Arrptr->Manningsn != NULL)
		n = Arrptr->Manningsn[source_index];
	else
		n = Parptr->FPn;
	ps_info->g_friction_squared_FP[ws_index] = g * n * n;

	if (Arrptr->SGCManningsn != NULL)
		n = Arrptr->SGCManningsn[source_index] * Arrptr->SGCManningsn[source_index];
	else
	{
		int group = Arrptr->SGCgroup[source_index];
		n = SGCptr->SGCn[group]; // note: already squared
	}
	ps_info->g_friction_squared_SG[ws_index] = g * n;

	SubGridCellInfo ws_cell;
	CopyToSubGridCellInfo(grid_cols_padded, ws_index, x, y, &ps_info->ws_cell, Parptr, Arrptr, SGCptr);
}

void BoundaryConditionToWaterSource(const int bc_index, const int ws_index, const int x, const int y,
	const int grid_cols_padded, const int grid_cols,
	const NUMERIC_TYPE g,
	WaterSource *bc_info,
	Pars *Parptr, Arrays *Arrptr,
	BoundCs * BCptr, SGCprams *SGCptr)
{
	bc_info->Ident[ws_index] = BCptr->BC_Ident[bc_index];
	bc_info->Val[ws_index] = BCptr->BC_Val[bc_index];
	bc_info->timeSeries[ws_index] = BCptr->BC_TimeSeries[bc_index];

	bc_info->Q_FP_old[ws_index] = C(0.0); //BCptr->BC_Q_FP_old[bc_index];
	bc_info->Q_SG_old[ws_index] = C(0.0); //BCptr->BC_Q_SG_old[bc_index];;

	// reading from non padded grids, so use grid_cols
	int source_index = y * grid_cols + x;

	NUMERIC_TYPE n;
	if (Arrptr->Manningsn != NULL)
		n = Arrptr->Manningsn[source_index];
	else
		n = Parptr->FPn;
	bc_info->g_friction_squared_FP[ws_index] = g * n * n;

	if (Arrptr->SGCManningsn != NULL)
		n = Arrptr->SGCManningsn[source_index] * Arrptr->SGCManningsn[source_index];
	else
	{
		int group = Arrptr->SGCgroup[source_index];
		n = SGCptr->SGCn[group]; // note: already squared
	}
	bc_info->g_friction_squared_SG[ws_index] = g * n;

	SubGridCellInfo ws_cell;
	CopyToSubGridCellInfo(grid_cols_padded, ws_index, x, y, &bc_info->ws_cell, Parptr, Arrptr, SGCptr);
}

/// convert point sources to list of rows
void InitPointSourceInfoStructure(PointSourceRowList * ps_layout,
	const int grid_rows, const int grid_cols_padded, const int grid_cols,
	const NUMERIC_TYPE g,
	Pars *Parptr, Arrays *Arrptr,
	BoundCs * BCptr, SGCprams *SGCptr)
{
	// count point sources per row
	int * ps_row_count_tmp = (int*)memory_allocate(grid_rows * sizeof(int));
	memset(ps_row_count_tmp, 0, grid_rows * sizeof(int));
	for (int i = 0; i < BCptr->numPS; i++)
	{
		int y = BCptr->ypi[i];
		ps_row_count_tmp[y]++;
	}

	int ps_count_row_max = 0;
	for (int y = 0; y < grid_rows; y++)
	{
		if (ps_row_count_tmp[y] > ps_count_row_max)
			ps_count_row_max = ps_row_count_tmp[y];
	}

	int row_cols_padded = ps_count_row_max + (GRID_ALIGN_WIDTH - (ps_count_row_max % GRID_ALIGN_WIDTH)) % GRID_ALIGN_WIDTH;
	ps_layout->row_cols_padded = row_cols_padded;
	ps_layout->ps_row_count = (int*)memory_allocate(grid_rows * sizeof(int));
	AllocateWaterSource(grid_rows * row_cols_padded, &(ps_layout->ps_info));

	/// create a list of point source indexes per row
	/// this will enable looping over rows to build the point source structure
	/// the reason to loop over rows while creating the point source structure is to ensure the memory is 
	/// allocated on the correct CPU's memory in a multi CPU machine (NUMA)
	int * ps_index_lookup_tmp = (int*)memory_allocate(grid_rows * row_cols_padded * sizeof(int));

	int * row_cursor = (int*)memory_allocate(grid_rows * sizeof(int));
	memset(row_cursor, 0, grid_rows * sizeof(int));

	for (int ps_index = 0; ps_index < BCptr->numPS; ps_index++)
	{
		int x = BCptr->xpi[ps_index];
		int y = BCptr->ypi[ps_index];

		int ps_new_index = (row_cols_padded * y) + row_cursor[y];
		row_cursor[y]++;
		ps_index_lookup_tmp[ps_new_index] = ps_index;
		printf("Point Source id:%d (%d,%d)\n", ps_index, x, y);
	}

	// omp loop over rows while allocating, ensure correct CPU first touch on multi CPU machine
#pragma omp parallel for default(shared) schedule(static)
	for (int y = 0; y < grid_rows; y++)
	{
		ps_layout->ps_row_count[y] = ps_row_count_tmp[y];

		for (int ps_i = 0; ps_i < ps_row_count_tmp[y]; ps_i++)
		{
			// write to index
			int ws_index = y * row_cols_padded + ps_i;
			// read from index
			int ps_index = ps_index_lookup_tmp[ws_index];

			PointSourceToWaterSource(ps_index, ws_index, grid_cols_padded, grid_cols,
				g,
				&ps_layout->ps_info, Parptr, Arrptr, BCptr, SGCptr);
		}
	}
	memory_free(&ps_index_lookup_tmp);
	memory_free(&row_cursor);
}


void InitBoundaryConditionStructure(BoundaryCondition * boundaryCondition,
	const int grid_rows, const int grid_cols_padded, const int grid_cols,
	const NUMERIC_TYPE g,
	Pars *Parptr, Arrays *Arrptr,
	BoundCs * BCptr, SGCprams *SGCptr)
{
	const int numBCs = 2 * grid_cols + 2 * grid_rows;

	AllocateWaterSource(numBCs, &boundaryCondition->bc_info);

	// initialise boundary
	{
		//printf("Init North boundary\n");
		int BCi, i, j, source_index;
		for (BCi = 0; BCi < grid_cols; BCi++)
		{
			// N(j=0) edge
			j = 0;
			i = BCi;
			//source_index = i;// + j * grid_cols;
			BoundaryConditionToWaterSource(BCi, BCi, i, j, grid_cols_padded, grid_cols, g, &boundaryCondition->bc_info, Parptr, Arrptr, BCptr, SGCptr);
		}
		//printf("Init East boundary\n");
		for (BCi = grid_cols; BCi < grid_cols + grid_rows; BCi++)
		{
			// E edge
			j = BCi - grid_cols;
			i = grid_cols - 1;
			//source_index = i + j * grid_cols;
			BoundaryConditionToWaterSource(BCi, BCi, i, j, grid_cols_padded, grid_cols, g, &boundaryCondition->bc_info, Parptr, Arrptr, BCptr, SGCptr);
		}
		//printf("Init South boundary\n");
		for (BCi = grid_cols + grid_rows; BCi < 2 * grid_cols + grid_rows; BCi++)
		{
			// S(j=ysz-1) edge
			j = grid_rows - 1;
			i = (2 * grid_cols + grid_rows) - BCi - 1;
			//source_index = i + j * grid_cols;
			BoundaryConditionToWaterSource(BCi, BCi, i, j, grid_cols_padded, grid_cols, g, &boundaryCondition->bc_info, Parptr, Arrptr, BCptr, SGCptr);
		}
		//printf("Init West boundary\n");
		for (BCi = 2 * grid_cols + grid_rows; BCi < numBCs; BCi++)
		{
			// W edge
			j = (2 * grid_cols + 2 * grid_rows) - BCi - 1;
			i = 0;
			//source_index = i + j * grid_cols;
			BoundaryConditionToWaterSource(BCi, BCi, i, j, grid_cols_padded, grid_cols, g, &boundaryCondition->bc_info, Parptr, Arrptr, BCptr, SGCptr);
		}
	}

}

void InitWeirStructure(EWeirType weirType,
	WeirLayout * weirs, const int grid_rows, const int grid_cols, const int grid_cols_padded,
	const int * flow_Qx_lookup_grid_tmp, const int * flow_Qy_lookup_grid_tmp,
	const NUMERIC_TYPE g,
	States * Statesptr, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr, const int verbose)
{

	weirs->weir_Qx_row_count = (int*)memory_allocate(grid_rows * sizeof(int));
	weirs->weir_Qy_row_count = (int*)memory_allocate(grid_rows * sizeof(int));

	if (Statesptr->weirs != ON)
	{
		memset(weirs->weir_Qx_row_count, 0, grid_rows * sizeof(int));
		memset(weirs->weir_Qy_row_count, 0, grid_rows * sizeof(int));
		return;
	}

	if (verbose == ON)	printf("Init Weir...\t");

	// count point sources per row
	int * weir_Qx_row_count_tmp = (int*)memory_allocate(grid_rows * sizeof(int));
	memset(weir_Qx_row_count_tmp, 0, grid_rows * sizeof(int));
	int * weir_Qy_row_count_tmp = (int*)memory_allocate(grid_rows * sizeof(int));
	memset(weir_Qy_row_count_tmp, 0, grid_rows * sizeof(int));

	weirs->weir_count = Arrptr->weir_count;
	// just keep the non padded list of weirs, not aligned with the rows and indexes.
	weirs->Weir_hc = Arrptr->Weir_hc;
	weirs->Weir_Cd = Arrptr->Weir_Cd;
	weirs->Weir_m = Arrptr->Weir_m;
	weirs->Weir_w = Arrptr->Weir_w;
	weirs->Weir_Fixdir = Arrptr->Weir_Fixdir;
	weirs->Weir_Typ = Arrptr->Weir_Typ;
	if (verbose == ON)
	{
		printf("\n");
	}
	// count qx weirs
	for (int y = 0; y < grid_rows; y++)
	{
		for (int x = 0; x < grid_cols - 1; x++)
		{
			int source_q_index = y * (grid_cols + 1) + x + 1;
			int weir_id = Arrptr->Weir_Identx[source_q_index];
			if (weir_id != -1 && weirType == weirs->Weir_Typ[weir_id])
			{
				weir_Qx_row_count_tmp[y]++;
				if (verbose == ON)
				{
					if (weirs->Weir_Typ[weir_id] == EWeir_Bridge)
						printf("Bridge [%d] (x) (%d,%d)\n", weir_id, x, y);
					else if (weirs->Weir_Typ[weir_id] == EWeir_Weir)
						printf("Weir [%d] (x) (%d,%d) %d\n", weir_id, x, y, weirs->Weir_Fixdir[weir_id]);
					else
						printf("Weir? [%d] (x) (%d,%d)\n", weir_id, x, y);
				}
			}
		}
	}

	//count qy weirs
	for (int y = 0; y < grid_rows - 1; y++)
	{
		for (int x = 0; x < grid_cols; x++)
		{
			int source_q_index = (y + 1) * (grid_cols + 1) + x;
			int weir_id = Arrptr->Weir_Identy[source_q_index];
			if (weir_id != -1 && weirType == weirs->Weir_Typ[weir_id])
			{
				weir_Qy_row_count_tmp[y]++;
				if (verbose == ON)
				{
					int weir_id = Arrptr->Weir_Identy[source_q_index];
					if (weirs->Weir_Typ[weir_id] == EWeir_Bridge)
						printf("Bridge [%d] (y) (%d,%d)\n", weir_id, x, y);
					else if (weirs->Weir_Typ[weir_id] == EWeir_Weir)
						printf("Weir [%d] (y) (%d,%d) %d\n", weir_id, x, y, weirs->Weir_Fixdir[weir_id]);
					else
						printf("Weir? [%d] (y) (%d,%d)\n", weir_id, x, y);
				}
			}
		}
	}

	int weir_count_row_max = 0;
	for (int y = 0; y < grid_rows; y++)
	{
		if (weir_Qx_row_count_tmp[y] > weir_count_row_max)
			weir_count_row_max = weir_Qx_row_count_tmp[y];
		if (weir_Qy_row_count_tmp[y] > weir_count_row_max)
			weir_count_row_max = weir_Qy_row_count_tmp[y];
	}

	int row_cols_padded = weir_count_row_max + (GRID_ALIGN_WIDTH - (weir_count_row_max % GRID_ALIGN_WIDTH)) % GRID_ALIGN_WIDTH;
	weirs->row_cols_padded = row_cols_padded;

	AllocateWeir(grid_rows * row_cols_padded, weirs);

	int total_weirs = 0;
	int total_bridges = 0;
	int total_other = 0;

	// omp loop over rows while allocating, ensure correct CPU first touch on multi CPU machine
#pragma omp parallel for default(shared) schedule(static)
	for (int y = 0; y < grid_rows; y++)
	{
		int row_weir_start = y * row_cols_padded;
		int weir_row_index;
		weirs->weir_Qx_row_count[y] = weir_Qx_row_count_tmp[y];
		weir_row_index = 0;
		for (int x = 0; x < grid_cols - 1; x++)
		{
			int source_index = y * grid_cols + x;
			int source_q_index = y * (grid_cols + 1) + x + 1;
			int weir_id = Arrptr->Weir_Identx[source_q_index];
			if (weir_id != -1 && weirType == weirs->Weir_Typ[weir_id])
			{
				weirs->weir_index_qx[row_weir_start + weir_row_index] = weir_id;
				weirs->Weir_grid_index[weir_id] = y * grid_cols_padded + x;

				weirs->Weir_Q_old_SG[weir_id] = Arrptr->QxSGold[source_q_index]; // copy from saved q if any

				int padded_q_index = ((y * grid_cols_padded) + x) + 1;

				weirs->Weir_pair_stream_flow_index[weir_id * 2] = flow_Qx_lookup_grid_tmp[padded_q_index - 1]; // flow to the west of the weir
				weirs->Weir_pair_stream_flow_index[weir_id * 2 + 1] = flow_Qx_lookup_grid_tmp[padded_q_index + 1]; // flow to the east of the weir

				if (weirs->Weir_Typ[weir_id] == EWeir_Bridge)
				{
					if (weirs->Weir_pair_stream_flow_index[weir_id * 2] == -1 ||
						weirs->Weir_pair_stream_flow_index[weir_id * 2 + 1] == -1)
					{
						// bridges use the upstream flow, which is assumed to be in the cell to the east or west
						printf("\nInvalid bridge cell. Bridge must have sub grid flows on either side.\n");
						exit(-1);
					}
					total_bridges++;
#ifdef _DEBUG
					//printf("Bridge (x) (%d,%d)\n", x, y);
#endif
				}
				else if (weirs->Weir_Typ[weir_id] == EWeir_Weir)
				{
#ifdef _DEBUG
					//printf("Weir (x) (%d,%d)\n", x, y);
#endif
					total_weirs++;
				}
				else
				{
					printf("Weir (x) Unknown (%d,%d)\n", x, y);
					total_other++;
				}

				int group0 = Arrptr->SGCgroup[source_index];
				int group1 = Arrptr->SGCgroup[source_index + 1];
				weirs->Weir_g_friction_sq[weir_id] = g * C(0.5)* (SGCptr->SGCn[group0] + SGCptr->SGCn[group1]);

				CopyToSubGridCellInfo(grid_cols_padded, 2 * weir_id, x, y, &weirs->cell_pair, Parptr, Arrptr, SGCptr);
				CopyToSubGridCellInfo(grid_cols_padded, 2 * weir_id + 1, x + 1, y, &weirs->cell_pair, Parptr, Arrptr, SGCptr);

				weir_row_index++;
			}
		}

		weirs->weir_Qy_row_count[y] = weir_Qy_row_count_tmp[y]; // last row is set to zero
		if (y < grid_rows - 1) // if not last row
		{
			weir_row_index = 0;
			for (int x = 0; x < grid_cols; x++)
			{
				int source_index = y * grid_cols + x;
				int source_q_index = (y + 1) * (grid_cols + 1) + x;
				int weir_id = Arrptr->Weir_Identy[source_q_index];
				if (weir_id != -1 && weirType == weirs->Weir_Typ[weir_id])
				{
					weirs->weir_index_qy[row_weir_start + weir_row_index] = weir_id;

					weirs->Weir_grid_index[weir_id] = y * grid_cols_padded + x;

					weirs->Weir_Q_old_SG[weir_id] = Arrptr->QySGold[source_q_index]; // copy from saved q if any

					int padded_q_index = ((y + 1) * grid_cols_padded) + x;

					weirs->Weir_pair_stream_flow_index[weir_id * 2] = flow_Qy_lookup_grid_tmp[padded_q_index - grid_cols_padded];
					weirs->Weir_pair_stream_flow_index[weir_id * 2 + 1] = flow_Qy_lookup_grid_tmp[padded_q_index + grid_cols_padded];

					if (weirs->Weir_Typ[weir_id] == EWeir_Bridge)
					{
						if (weirs->Weir_pair_stream_flow_index[weir_id * 2] == -1 ||
							weirs->Weir_pair_stream_flow_index[weir_id * 2 + 1] == -1)
						{
							// bridges use the upstream flow, which is assumed to be in the cell to the north or south
							printf("\nInvalid bridge cell. Bridge must have sub grid flows on either side.\n");
							exit(-1);
						}
						total_bridges++;
#ifdef _DEBUG
						//printf("Bridge (y) (%d,%d)\n", x, y);
#endif
					}
					else if (weirs->Weir_Typ[weir_id] == EWeir_Weir)
					{
#ifdef _DEBUG
						//printf("Weir (y) (%d,%d)\n", x, y);
#endif
						total_weirs++;
					}
					else
					{
						printf("Weir (y) Unknown (%d,%d)\n", x, y);
						total_other++;
					}

					int group0 = Arrptr->SGCgroup[source_index];
					int group1 = Arrptr->SGCgroup[source_index + grid_cols];
					weirs->Weir_g_friction_sq[weir_id] = g*C(0.5)* (SGCptr->SGCn[group0] + SGCptr->SGCn[group1]);

					CopyToSubGridCellInfo(grid_cols_padded, 2 * weir_id, x, y, &weirs->cell_pair, Parptr, Arrptr, SGCptr);
					CopyToSubGridCellInfo(grid_cols_padded, 2 * weir_id + 1, x, y + 1, &weirs->cell_pair, Parptr, Arrptr, SGCptr);

					weir_row_index++;
				}
			}
		}
	}

	if (verbose == ON)
	{
		if (total_weirs > 0)
			printf("Weirs: %d\n", total_weirs);
		if (total_bridges > 0)
			printf("Bridges: %d\n", total_bridges);
	}
}

void InitRoutingStructure_row(NUMERIC_TYPE * route_V_ratio_per_sec_qx, NUMERIC_TYPE * route_V_ratio_per_sec_qy,
	const int j, const int not_last_row,
	const int grid_rows, const int grid_cols, const int grid_cols_padded,
	Arrays *Arrptr)
{
	memset(route_V_ratio_per_sec_qx + (j * grid_cols_padded), 0, sizeof(NUMERIC_TYPE) * grid_cols_padded);
	memset(route_V_ratio_per_sec_qy + ((j + 1) * grid_cols_padded), 0, sizeof(NUMERIC_TYPE) * grid_cols_padded);
	for (int i = 0; i < grid_cols - 1; i++)
	{
		int source_index = i + j * grid_cols; //read from old non-padded data
		int source_index_next = source_index + 1; //read from old non-padded data
		int dest_q_index = (i + j * grid_cols_padded) + 1; // write to padded data

		if (Arrptr->DEM[source_index] > Arrptr->DEM[source_index_next] &&
			Arrptr->FlowDir[source_index] == source_index_next)
		{
			route_V_ratio_per_sec_qx[dest_q_index] = -1 / Arrptr->RouteInt[source_index];
		}
		else if (Arrptr->DEM[source_index] < Arrptr->DEM[source_index_next] &&
			Arrptr->FlowDir[source_index_next] == source_index)
		{
			route_V_ratio_per_sec_qx[dest_q_index] = 1 / Arrptr->RouteInt[source_index_next];
		}
	}
	if (not_last_row)
	{
		for (int i = 0; i < grid_cols; i++)
		{
			int source_index = i + j * grid_cols; //read from old non-padded data
			int source_index_next = source_index + grid_cols; //read from old non-padded data
			int dest_q_index = (i + j * grid_cols_padded) + grid_cols_padded; // write to padded data

			if (Arrptr->DEM[source_index] > Arrptr->DEM[source_index_next] &&
				Arrptr->FlowDir[source_index] == source_index_next)
			{
				route_V_ratio_per_sec_qy[dest_q_index] = -1 / Arrptr->RouteInt[source_index];
			}
			else if (Arrptr->DEM[source_index] < Arrptr->DEM[source_index_next] &&
				Arrptr->FlowDir[source_index_next] == source_index)
			{
				route_V_ratio_per_sec_qy[dest_q_index] = 1 / Arrptr->RouteInt[source_index_next];
			}
		}
	}
}

void InitDamStructure(const int grid_rows, const int grid_cols, const int grid_cols_padded,	States *Statesptr, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr, DamData *Damptr, int verbose)
{
	int n = 0;
	int nt = 0;
	int gr;
	Damptr->TotalEdge = C(0.0);
	for (int d = 0; d < Damptr->NumDams; d++)
	{
		Damptr->TotalEdge += Damptr->Edgenos[d];
		Damptr->DamVol[d] = (Damptr->Volmax[d] / Damptr->DamHeight[d])*(Damptr->InitialHeight[d] - (Damptr->SpillHeight[d] - Damptr->DamHeight[d]));
	}
	Damptr->DynamicEdge = new int [3*Damptr->TotalEdge]();
//	Damptr->TotalEdge = nt;
	Damptr->DynamicEdgeData = memory_allocate_zero_numeric_legacy(12*Damptr->TotalEdge);
	for (int i = 0; i < grid_cols; i++) for (int j = 0; j < grid_rows; j++) 
	{
		if (Arrptr->DamMask[i + j*grid_cols] > C(0.0))
		{
			int source_index = i + j*grid_cols;
			NUMERIC_TYPE cell_lenght = (Arrptr->dx[j*grid_cols] + Arrptr->dy[j*grid_cols]) / C(2.0);
			//cell_lenght = cell_lenght / C(2.0);
			gr = Arrptr->SGCgroup[source_index];

			Damptr->DynamicEdge[2 * n] = i + j*grid_cols_padded;
			Damptr->DynamicEdge[(2 * n) + 1] = Arrptr->DamMask[i + j*grid_cols];
			Damptr->DynamicEdge[(2 * n) + 2] = gr;
			Damptr->DynamicEdgeData[12 * n] = 0; //Qold
			Damptr->DynamicEdgeData[12 * n + 1] = 0; // Q combined
			Damptr->DynamicEdgeData[12 * n + 2] = 0; // Wfp
			Damptr->DynamicEdgeData[12 * n + 3] = Arrptr->Manningsn[source_index]; //nfp
			Damptr->DynamicEdgeData[12 * n + 4] = Arrptr->DEM[source_index]; // Zfp
			Damptr->DynamicEdgeData[12 * n + 5] = SGCptr->SGCn[gr]; // nSGC
			Damptr->DynamicEdgeData[12 * n + 6] = Arrptr->SGCz[source_index]; // BankFullDepth
			Damptr->DynamicEdgeData[12 * n + 7] = 0; // QoldSGC
			Damptr->DynamicEdgeData[12 * n + 8] = cell_lenght; // cell_lenght
			Damptr->DynamicEdgeData[12 * n + 9] = SGCptr->SGCs[gr]; // sSGC
			Damptr->DynamicEdgeData[12 * n + 10] = SGCptr->SGCchantype[gr]; // ChantypeSGC
			Damptr->DynamicEdgeData[12 * n + 11] = Arrptr->SGCwidth[source_index];; // widthSGC

//			printf("%d\n", n);
			if (i > C(0.0) && Arrptr->DamMask[i-1 + j*grid_cols] < C(0.0))
			{
				Damptr->DynamicEdgeData[12 * n + 2] += Arrptr->dx[j*grid_cols];
			}
			if (i < grid_cols && Arrptr->DamMask[i + 1 + j*grid_cols] < C(0.0))
			{
				Damptr->DynamicEdgeData[12 * n + 2] += Arrptr->dx[j*grid_cols];
			}
			if (j > C(0.0) && Arrptr->DamMask[i + (j - 1)*grid_cols] < C(0.0))
			{
				Damptr->DynamicEdgeData[12 * n + 2] += Arrptr->dy[(j-1)*grid_cols];
			}
			if (j < grid_rows && Arrptr->DamMask[i + (j + 1)*grid_cols] < C(0.0))
			{
				Damptr->DynamicEdgeData[12 * n + 2] += Arrptr->dy[(j+1)*grid_cols];
			}		
			n++;
		}
	}
}


void Fast_MainInit(Fnames *Fnameptr, Files *Fptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr, Stage *Locptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, SGCprams *SGCptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, DamData *Damptr, int verbose)
{
	const int grid_cols = Parptr->xsz;
	const int grid_rows = Parptr->ysz;
	// keep consistent stride throughout
	// extra 1 column needed for qx, and friction
	// extra padding every row has at least 64 bytes of blank padding on the right
	int grid_cols_padded = grid_cols + 1 + (64 / sizeof(NUMERIC_TYPE));
	grid_cols_padded += (GRID_ALIGN_WIDTH - (grid_cols_padded % GRID_ALIGN_WIDTH)) % GRID_ALIGN_WIDTH;

	// clean up memory that is not required

#ifndef RESULT_CHECK
	//memory_free_legacy(&Arrptr->DEM ); 
	//memory_free_legacy(&Arrptr->H );
	memory_free_legacy(&Arrptr->Qx); //todo copy in case of checkpoint
	memory_free_legacy(&Arrptr->Qy); //todo copy in case of checkpoint
	memory_free_legacy(&Arrptr->Qxold); //todo copy in case of checkpoint?
	memory_free_legacy(&Arrptr->Qyold); //todo copy in case of checkpoint?
	memory_free_legacy(&Arrptr->U);
	//memory_free_legacy(&Arrptr->V);

	/* TRENT additions  */
	memory_free_legacy(&Arrptr->HU);
	memory_free_legacy(&Arrptr->HV);
	memory_free_legacy(&Arrptr->RSHU);
	memory_free_legacy(&Arrptr->LSHU);
	memory_free_legacy(&Arrptr->RSHV);
	memory_free_legacy(&Arrptr->LSHV);
	memory_free_legacy(&Arrptr->BSHU);
	memory_free_legacy(&Arrptr->TSHU);
	memory_free_legacy(&Arrptr->BSHV);
	memory_free_legacy(&Arrptr->TSHV);
	memory_free_legacy(&Arrptr->FHx);
	memory_free_legacy(&Arrptr->FHUx);
	memory_free_legacy(&Arrptr->FHVx);
	memory_free_legacy(&Arrptr->FHy);
	memory_free_legacy(&Arrptr->FHUy);
	memory_free_legacy(&Arrptr->FHVy);

	//memory_free_legacy(&Arrptr->FlowDir);
	memory_free_legacy(&Arrptr->Route_dH);
	//memory_free_legacy(&Arrptr->RouteInt);

	/* ---------------- */
	memory_free_legacy(&Arrptr->maxH); //todo copy in case of checkpoint
	memory_free_legacy(&Arrptr->maxHtm); //todo copy in case of checkpoint
	memory_free_legacy(&Arrptr->initHtm); //todo copy in case of checkpoint
	memory_free_legacy(&Arrptr->totalHtm); //todo copy in case of checkpoint
	//memory_free_legacy(&Arrptr->Manningsn );
	//memory_free_legacy(&Arrptr->SGCManningsn );
	memory_free_legacy(&Arrptr->paerial);
	memory_free_legacy(&Arrptr->pbound);
	//memory_free_legacy(&Arrptr->Weir_hc); // lists of weirs used
	//memory_free_legacy(&Arrptr->Weir_Cd);
	//memory_free_legacy(&Arrptr->Weir_m);
	//memory_free_legacy(&Arrptr->Weir_w);
	//memory_free_legacy(&Arrptr->Weir_Fixdir);
	//memory_free_legacy(&Arrptr->Weir_Typ);

	//memory_free_legacy(&Arrptr->evap );
	//memory_free_legacy(&Arrptr->rain );

	//memory_free_legacy(&Arrptr->ChanMask);
	memory_free_legacy(&Arrptr->SegMask);

	memory_free_legacy(&Arrptr->TRecx);
	memory_free_legacy(&Arrptr->TRecy);

	memory_free_legacy(&Arrptr->LimQx);
	memory_free_legacy(&Arrptr->LimQy);
	memory_free_legacy(&Arrptr->Vx);
	memory_free_legacy(&Arrptr->Vy);
	memory_free_legacy(&Arrptr->maxVx);
	memory_free_legacy(&Arrptr->maxVy);
	memory_free_legacy(&Arrptr->Vc);
	memory_free_legacy(&Arrptr->maxVc);
	memory_free_legacy(&Arrptr->maxVcH);
	memory_free_legacy(&Arrptr->maxHaz);
	//memory_free_legacy(&Arrptr->SGCwidth );
	//memory_free_legacy(&Arrptr->SGCz);

	//memory_free_legacy(&Arrptr->SGCbfH );
	//memory_free_legacy(&Arrptr->SGCVol );
	memory_free_legacy(&Arrptr->SGCdVol);
	//memory_free_legacy(&Arrptr->SGCbfV );
	//memory_free_legacy(&Arrptr->SGCc );
	memory_free_legacy(&Arrptr->SGCFlowWidth);
	memory_free_legacy(&Arrptr->SGCdx);
	memory_free_legacy(&Arrptr->SGCcat_area);
	//memory_free_legacy(&Arrptr->dx );
	//memory_free_legacy(&Arrptr->dy );
	//memory_free_legacy(&Arrptr->dA ); 

	//memory_free_legacy(&Arrptr->SGCgroup);

#endif

	const NUMERIC_TYPE g = Solverptr->g;

	TimeSeries * evap_time_series = Arrptr->evap;
	TimeSeries * rain_time_series = Arrptr->rain;

	NUMERIC_TYPE * rain_grid       = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	NUMERIC_TYPE * dist_infil_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows); // grid for dirtributed infiltration data

	NetCDFVariable evap_grid;
#if (_NETCDF == 1)
	if (Statesptr->calc_evap == TIME_SPACE) {
		read_file_netCDF_start(Fnameptr->evapfilename, "evaporation", &evap_grid);
	}
#endif
	evap_grid.data = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);

	// digital elevation model (z)
	NUMERIC_TYPE * dem_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	
	//dx_col: read only data - column vector - one cell stores the grid dx for every cell in the row
	// cell x may vary with latitude on large models
	NUMERIC_TYPE *dx_col = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * (grid_rows + 1));

	//dy_col: read only data - column vector - one cell stores the grid dx for every cell in the row
	// cell y may vary with latitude on large models
	NUMERIC_TYPE *dy_col = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * (grid_rows + 1));

	//dA: read only data - column vector - one cell stores the grid dx for every cell in the row
	// cell area may vary with latitude on large models
	NUMERIC_TYPE *cell_area_col = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * (grid_rows + 1));

	/// friction between the indexed cell and the next cell (right)
	// row 0 has friction between row 0 and row 1
	// data is offset to the right by 1 column for vectorization alignment with qx
	NUMERIC_TYPE * g_friction_sq_x_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	NUMERIC_TYPE * friction_x_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);

	/// friction between the indexed cell and the next cell (below)
	/// column 0 has friction between column 0 and column 1
	// data is offset to the down by 1 row for vectorization alignment with qy
	NUMERIC_TYPE * g_friction_sq_y_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * (grid_rows + 1));
	NUMERIC_TYPE * friction_y_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * (grid_rows + 1));
	
	NUMERIC_TYPE * h_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	NUMERIC_TYPE * volume_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	//int * ChanMask_grid = (int*)memory_allocate(sizeof(int) * grid_cols_padded * grid_rows);

	NUMERIC_TYPE *maxH_grid, *maxHtm_grid, *initHtm_grid, *totalHtm_grid, *maxVc_grid, *maxVc_height_grid, *maxHazard_grid;
	maxH_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	maxHtm_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	initHtm_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	totalHtm_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	maxVc_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	maxVc_height_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	maxHazard_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);

	// qx is in m3/s
	// column 0 is flow to/from out of bounds into the domain
	// last column is flow to/from out of bounds into the domain
	// this array has one more column than the h/dem etc
	NUMERIC_TYPE* Qx_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	NUMERIC_TYPE* Qx_old_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	NUMERIC_TYPE* Vx_max_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	NUMERIC_TYPE* Vx_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	// floodplain ywidth is defined on Qx grid
	NUMERIC_TYPE* Fp_ywidth = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);

	// qy is in m3/s
	// row 0 is flow to/from out of bounds into the domain
	// the last row is flow to/from out of bounds into the domain
	//note extra row on Qy
	NUMERIC_TYPE* Qy_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * (grid_cols_padded * (grid_rows + 1)));
	NUMERIC_TYPE* Qy_old_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * (grid_cols_padded * (grid_rows + 1)));
	NUMERIC_TYPE* Vy_max_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * (grid_cols_padded * (grid_rows + 1)));
	NUMERIC_TYPE* Vy_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * (grid_cols_padded * (grid_rows + 1)));
	// floodplain xwidth is defined on Qy grid
	NUMERIC_TYPE* Fp_xwidth = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * (grid_rows + 1));

	// sub-grid
	//NUMERIC_TYPE *SGC_width_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	NUMERIC_TYPE *SGC_BankFullHeight_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	//NUMERIC_TYPE *SGC_BankFullVolume_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	//NUMERIC_TYPE *SGC_c_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	//NUMERIC_TYPE *SGC_Qx_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	//NUMERIC_TYPE *SGC_Qy_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * (grid_rows + 1));
	//flood plain flow into cell needs to be removed for the region above the channel (Neil 2012 Figure 1 (C) )
	//NUMERIC_TYPE *SGC_ChannelWidth_X_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	//flood plain flow into cell needs to be removed for the region above the channel (Neil 2012 Figure 1 (C) )
	//NUMERIC_TYPE *SGC_ChannelWidth_Y_grid = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * (grid_rows + 1));

	//NUMERIC_TYPE * SGC_mannings;
	//if (Arrptr->SGCManningsn != NULL)
	//	SGC_mannings = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded * grid_rows);
	//else
	//	SGC_mannings = NULL;

	// count the number of threads
	// omp_get_num_threads() - current number of threads can't be used outside of pragma omp
	// omp_get_max_threads() - is not consistent across implementations
	// most reliable way is to just count them
	int thread_count = 0;
#pragma omp parallel reduction (+:thread_count)
	{
		thread_count++;
	}
	if (verbose == ON) printf("OMP thread count: %d\n", thread_count);

#ifdef __unix__
	int* numa_nodes = (int*)memory_allocate(sizeof(int)*thread_count);
#pragma omp parallel default(shared)
	{
		int node = numa_node_of_cpu(sched_getcpu());

		numa_nodes[omp_get_thread_num()] = node;
	}
	printf("NUMA:\n");
	for (int i = 0; i < thread_count; i++)
	{
		printf("|%d", numa_nodes[i]);
	}
	printf("|\n");
	memory_free(&numa_nodes);
#endif

	int block_count = thread_count;

	WetDryRowBound* wet_dry_bounds = new WetDryRowBound();
	AllocateWetDryRowBound(grid_rows, block_count, wet_dry_bounds);

	for (int j = 0; j < grid_rows; j++)
	{
		wet_dry_bounds->fp_h[j].start = 0;
		wet_dry_bounds->fp_h[j].end = grid_cols;
	}

	//int *SGC_group_grid = (int*)memory_allocate(sizeof(int) * grid_cols_padded * grid_rows);
#ifdef _DEBUG
	printf("InitPointSourceInfoStructure\n");
#endif
	PointSourceRowList * ps_layout = new PointSourceRowList();
	InitPointSourceInfoStructure(ps_layout, grid_rows, grid_cols_padded, grid_cols, g, Parptr, Arrptr, BCptr, SGCptr);
#ifdef _DEBUG
	printf("InitBoundaryConditionStructure\n");
#endif
	BoundaryCondition * boundaryCondition = new BoundaryCondition();
	InitBoundaryConditionStructure(boundaryCondition, grid_rows, grid_cols_padded, grid_cols, g, Parptr, Arrptr, BCptr, SGCptr);

	SubGridRowList * sub_grid_layout_rows = new SubGridRowList();
	SubGridState * sub_grid_state_rows = new SubGridState();
	SubGridRowList * sub_grid_layout_blocks = new SubGridRowList();
	SubGridState * sub_grid_state_blocks = new SubGridState();

	//temporary data for building weir (bridges) structure
	int* flow_Qx_lookup_grid_tmp = (int*)memory_allocate(sizeof(int)* (grid_rows + 1) * grid_cols_padded);
	SetArrayValue(flow_Qx_lookup_grid_tmp, -1, (grid_rows + 1) * grid_cols_padded);
	//temporary data for building weir (bridges) structure
	int* flow_Qy_lookup_grid_tmp = (int*)memory_allocate(sizeof(int)* (grid_rows + 1) * grid_cols_padded);
	SetArrayValue(flow_Qy_lookup_grid_tmp, -1, (grid_rows + 1) * grid_cols_padded);

	// PFU initialise floodplain width variables (needed before setting subgrid)


#ifdef _DEBUG
	printf("InitSubGridStructure\n");
#endif
	InitSubGridStructureByRows(sub_grid_state_rows, sub_grid_layout_rows,
		grid_cols_padded,
		g,
		flow_Qx_lookup_grid_tmp, flow_Qy_lookup_grid_tmp,
		Fp_xwidth,	Fp_ywidth,
		Statesptr, Parptr, Arrptr, SGCptr, verbose);
#if _SGM_BY_BLOCKS == 1
	InitSubGridStructureByBlocks(sub_grid_state_blocks, sub_grid_layout_blocks,
		block_count,
		grid_cols_padded,
		g,
		flow_Qx_lookup_grid_tmp, flow_Qy_lookup_grid_tmp,
		FP_xwidth, Fp_ywidth,
		Statesptr, Parptr, Arrptr, SGCptr, verbose);
#endif
#ifdef _DEBUG
	printf("InitWeirStructure\n");
#endif
	WeirLayout * weirs_weirs = new WeirLayout();
	WeirLayout * weirs_bridges = new WeirLayout();
	InitWeirStructure(EWeir_Weir, weirs_weirs, grid_rows, grid_cols, grid_cols_padded,
		flow_Qx_lookup_grid_tmp, flow_Qy_lookup_grid_tmp, g,
		Statesptr, Parptr, Arrptr, SGCptr, verbose);
	InitWeirStructure(EWeir_Bridge, weirs_bridges, grid_rows, grid_cols, grid_cols_padded,
		flow_Qx_lookup_grid_tmp, flow_Qy_lookup_grid_tmp, g,
		Statesptr, Parptr, Arrptr, SGCptr, verbose);
	
	memory_free(&flow_Qx_lookup_grid_tmp);
	memory_free(&flow_Qy_lookup_grid_tmp);

	// Added by FEOL for Dams
	if (Statesptr->DamMode == ON) InitDamStructure(grid_rows, grid_cols, grid_cols_padded, Statesptr, Parptr, Arrptr,SGCptr, Damptr,verbose);
	
	// Added by JCN for super grid channels
	SuperGridLinksList * Super_linksptr = new SuperGridLinksList();
	if (Statesptr->ChanMaskRead == ON)
	{ 
		InitSuperLinksStructure(grid_rows, grid_cols, grid_cols_padded, Super_linksptr, Statesptr, Parptr, Arrptr, SGCptr, Solverptr, Fnameptr, verbose);
		//memory_free_legacy(&Arrptr->ChanMask);
	}
	memory_free_legacy(&Arrptr->SGCz);

	RouteDynamicList * route_dynamic_list = new RouteDynamicList();
	AllocateRoutingDynamicList(grid_rows, grid_cols_padded, route_dynamic_list);

	NUMERIC_TYPE * route_V_ratio_per_sec_qx = NULL;
	NUMERIC_TYPE * route_V_ratio_per_sec_qy = NULL;
	if (Statesptr->routing == ON)
	{
		route_V_ratio_per_sec_qx = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_rows * grid_cols_padded);
		route_V_ratio_per_sec_qy = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * (grid_rows + 1) * grid_cols_padded);
	}

	/// set up temp data for use by each thread
	NUMERIC_TYPE ** tmp_thread_data = (NUMERIC_TYPE**)memory_allocate(sizeof(NUMERIC_TYPE*) * thread_count);
	for (int thread_id = 0; thread_id < thread_count; thread_id++)
	{
		//allocate 1 row of temp per thread
		tmp_thread_data[thread_id] = (NUMERIC_TYPE*)memory_allocate(sizeof(NUMERIC_TYPE) * grid_cols_padded);
	}

	//first touch and zero tmp_thread_data in the thread it will be used by
#pragma omp parallel default(shared)
	{
		int thread_id = omp_get_thread_num();
		memset(tmp_thread_data[thread_id], 0, sizeof(NUMERIC_TYPE) * grid_cols_padded);
	}


	//start with even load balance
	SGC2_UpdateLoadBalance(grid_rows, grid_cols_padded, sub_grid_layout_rows, wet_dry_bounds);

#pragma omp parallel for default(shared) schedule(static)
	//for (int j = 0; j < grid_rows; j++)
	for (int block_index = 0; block_index < block_count; block_index++)
	{
		const int start_y = wet_dry_bounds->block_row_bounds[block_index].start;
		const int end_y = wet_dry_bounds->block_row_bounds[block_index].end;

		for (int j = start_y; j < end_y; j++)
		{
#ifdef _DEBUG
			//printf("j %d thread %d\n", j, omp_get_thread_num());
#endif

			int source_row_index = j * grid_cols;
			int dest_row_index = j * grid_cols_padded;

			int padding_count = grid_cols_padded - grid_cols;

			int source_bytes_per_row = sizeof(NUMERIC_TYPE) * grid_cols;
			int dest_bytes_per_row = sizeof(NUMERIC_TYPE) * grid_cols_padded;
			int padding = sizeof(NUMERIC_TYPE) * padding_count;
			
			int source_bytes_per_row_int = sizeof(int)*grid_cols;
			int dest_bytes_per_row_int = sizeof(int)*grid_cols_padded;
			int padding_int = dest_bytes_per_row_int - source_bytes_per_row_int;

			int start = -1, end = -1;
			for (int i = 0; i < grid_cols; i++)
			{
				if (start == -1 && Arrptr->DEM[source_row_index + i] != DEM_NO_DATA)
					start = i;
				if (Arrptr->DEM[source_row_index + i] != DEM_NO_DATA)
					end = i;
			}
			wet_dry_bounds->dem_data[j].start = start;
			wet_dry_bounds->dem_data[j].end = min(grid_cols, end + 1);

			//initialize to bounds of DEM - some models do not cover the entire grid.
			wet_dry_bounds->fp_h[j].start = start;
			wet_dry_bounds->fp_h[j].end = wet_dry_bounds->dem_data[j].end;
			wet_dry_bounds->fp_h_prev[j].start = start;
			wet_dry_bounds->fp_h_prev[j].end = wet_dry_bounds->dem_data[j].end;
			wet_dry_bounds->fp_vol[j].start = start;
			wet_dry_bounds->fp_vol[j].end = wet_dry_bounds->dem_data[j].end;

			route_dynamic_list->row_route_qx_count[j] = 0;
			route_dynamic_list->row_route_qy_count[j] = 0;
			memset(route_dynamic_list->route_list_i_lookup_qx + dest_row_index, 0, sizeof(int) * grid_cols_padded);
			memset(route_dynamic_list->route_list_i_lookup_qy + dest_row_index, 0, sizeof(int) * grid_cols_padded);

			memcpy(h_grid + dest_row_index, Arrptr->H + source_row_index, source_bytes_per_row);
			memset(h_grid + dest_row_index + grid_cols, 0, padding);

			// todo: check initial volume is consistent with height - initial code reading suggests that it is
			memcpy(volume_grid + dest_row_index, Arrptr->SGCVol + source_row_index, source_bytes_per_row);
			memset(volume_grid + dest_row_index + grid_cols, 0, padding);

			memcpy(dem_grid + dest_row_index, Arrptr->DEM + source_row_index, source_bytes_per_row);
			memset(dem_grid + dest_row_index + grid_cols, 0, padding);

			//Initiate distributed rainfall mask ASmith 2015
			if (Statesptr->rainfallmask == ON)
			{
				memcpy(rain_grid + dest_row_index, Arrptr->Rainmask + source_row_index, source_bytes_per_row);
				memset(rain_grid + dest_row_index + grid_cols, 0, padding);
			}
			if (Statesptr->calc_distributed_infiltration == ON)
			{
				memcpy(dist_infil_grid + dest_row_index, Arrptr->dist_infiltration+ source_row_index, source_bytes_per_row);
				memset(dist_infil_grid + dest_row_index + grid_cols, 0, padding);
			}

			//Initiate distributed evaporation grid
			if (Statesptr->calc_evap == TIME_SPACE)
			{
				memset(evap_grid.data + dest_row_index, 0, sizeof(NUMERIC_TYPE) * grid_cols_padded);
			}

			/// sub grid - grid version
			//memcpy(SGC_width_grid + dest_row_index, Arrptr->SGCwidth + source_row_index, source_bytes_per_row);
			//memset(SGC_width_grid + dest_row_index + grid_cols, 0, padding);
			memcpy(SGC_BankFullHeight_grid + dest_row_index, Arrptr->SGCbfH + source_row_index, source_bytes_per_row);
			memset(SGC_BankFullHeight_grid + dest_row_index + grid_cols, 0, padding);

			//memcpy(ChanMask_grid + dest_row_index, Arrptr->ChanMask + source_row_index, sizeof(int) * grid_cols);
			//memset(ChanMask_grid + dest_row_index + grid_cols, 0, sizeof(int) * padding_count);

			//memcpy(SGC_BankFullVolume_grid + dest_row_index, Arrptr->SGCbfV + source_row_index, source_bytes_per_row);
			//memset(SGC_BankFullVolume_grid + dest_row_index + grid_cols, 0, padding);

			//memcpy(SGC_c_grid + dest_row_index, Arrptr->SGCc + source_row_index, source_bytes_per_row);
			//memset(SGC_c_grid + dest_row_index + grid_cols, 0, padding);

			//basic mannings copy - improve with 1d sub-grid-channel
			//if (SGC_mannings != NULL)
			//{
			//	memcpy(SGC_mannings + dest_row_index, Arrptr->SGCManningsn + source_row_index, source_bytes_per_row);
			//	memset(SGC_mannings + dest_row_index + grid_cols, 0, padding);
			//}

			//memcpy(SGC_group_grid + dest_row_index, Arrptr->SGCgroup + source_row_index, sizeof(int) * grid_cols);
			//memset(SGC_group_grid + dest_row_index + grid_cols, 0, sizeof(int) * padding_count);
			// end sub grid		

			//copy only the first cooumn of the cell width and area data.
			//data is constant across the x dimention
			dx_col[j] = Arrptr->dx[source_row_index];
			dy_col[j] = Arrptr->dy[source_row_index];
			cell_area_col[j] = Arrptr->dA[source_row_index];
			//memcpy(dx_col + dest_row_index, Arrptr->dx + source_row_index, source_bytes_per_row);
			//memcpy(dy_col + dest_row_index, Arrptr->dy + source_row_index, source_bytes_per_row);

			bool not_last_row = (j < grid_rows - 1);

			//memset(maxH_grid + dest_row_index, 0, dest_bytes_per_row);
			memset(totalHtm_grid + dest_row_index, 0, dest_bytes_per_row);
			memset(maxVc_grid + dest_row_index, 0, dest_bytes_per_row);
			memset(maxVc_height_grid + dest_row_index, 0, dest_bytes_per_row);
			memset(maxHazard_grid + dest_row_index, 0, dest_bytes_per_row);

			// ensure first value is zero
			g_friction_sq_x_grid[dest_row_index] = 0;
			friction_x_grid[dest_row_index] = 0;
			for (int i = 0; i < grid_cols; i++)
			{
				int dest_index = dest_row_index + i;
				int source_index = source_row_index + i;
				// update the h value to be always relative to the dem z
				// in the previous version h was conditionally h above dem z or h above SGCz
				if (SGC_BankFullHeight_grid[dest_index] == NULLVAL)
					SGC_BankFullHeight_grid[dest_index] = C(0.0);
				else
					h_grid[dest_index] -= SGC_BankFullHeight_grid[dest_index];

				NUMERIC_TYPE fn_x, fn_y;
				if (Arrptr->Manningsn != NULL)
				{
					fn_x = C(0.5)*(Arrptr->Manningsn[source_index] + Arrptr->Manningsn[source_index + 1]);
				}
				else
				{
					fn_x = Parptr->FPn;
				}
				//offset by one column to match qx alignment
				friction_x_grid[dest_index + 1] = fn_x;
				g_friction_sq_x_grid[dest_index + 1] = g * fn_x * fn_x;

				if (not_last_row)
				{
					if (Arrptr->Manningsn != NULL)
					{
						fn_y = C(0.5)*(Arrptr->Manningsn[source_index] + Arrptr->Manningsn[source_index + grid_cols]);
					}
					else
					{
						fn_y = Parptr->FPn;
					}
					//offset by one row to match qx alignment
					friction_y_grid[dest_index + grid_cols_padded] = fn_y;
					g_friction_sq_y_grid[dest_index + grid_cols_padded] = g * fn_y * fn_y;
				}
				else
				{
					//store an additional cell dx and dy for the additional row of flow
					dx_col[j + 1] = Arrptr->dx[source_row_index];
					dy_col[j + 1] = Arrptr->dy[source_row_index];
					cell_area_col[j + 1] = Arrptr->dA[source_row_index];
				}

				maxHtm_grid[dest_index] = (NUMERIC_TYPE)NULLVAL;
				initHtm_grid[dest_index] = (NUMERIC_TYPE)NULLVAL;
				maxH_grid[dest_index] = C(0.0);// -SGC_BankFullHeight_grid[dest_index];
			}
			for (int i = grid_cols; i < grid_cols_padded; i++)
			{
				int dest_index = dest_row_index + i;
				maxH_grid[dest_index] = 0;
				maxHtm_grid[dest_index] = 0;
				initHtm_grid[dest_index] = 0;
				g_friction_sq_x_grid[dest_index] = 0;
				g_friction_sq_y_grid[grid_cols_padded + dest_index] = 0;
				friction_x_grid[dest_index] = 0;
				friction_y_grid[grid_cols_padded + dest_index] = 0;
			}

			// TODO q copy from source (used as oldq) - for checkpoint load
			memset(Qx_old_grid + dest_row_index, 0, dest_bytes_per_row);
			memset(Qx_grid + dest_row_index, 0, dest_bytes_per_row);
			memset(Vx_max_grid + dest_row_index, 0, dest_bytes_per_row);
			memset(Vx_grid + dest_row_index, 0, dest_bytes_per_row);
			//memset(Fp_xwidth + dest_row_index, 0, dest_bytes_per_row); PFU needs to be initialised before init subgrid

			if (j == 0)
			{
				for (int i = 0; i < grid_cols_padded; i++)
				{
					g_friction_sq_y_grid[i] = 0;
				}
				// TODO q copy from source (used as oldq) - for checkpoint load
				memset(Qy_old_grid, 0, dest_bytes_per_row);
				memset(Qy_grid, 0, dest_bytes_per_row);
				memset(Vy_max_grid, 0, dest_bytes_per_row);
				memset(Vy_grid, 0, dest_bytes_per_row);
				memset(g_friction_sq_y_grid, 0, dest_bytes_per_row);
				memset(friction_y_grid, 0, dest_bytes_per_row);
				//memset(Fp_ywidth, 0, dest_bytes_per_row); PFU needs to be initialised before init subgrid

				if (Statesptr->routing == ON)
					memset(route_V_ratio_per_sec_qy, 0, dest_bytes_per_row);
			}
			//next row
			// TODO q copy from source (used as oldq) - for checkpoint load
			memset(Qy_old_grid + dest_row_index + grid_cols_padded, 0, dest_bytes_per_row);
			memset(Qy_grid + dest_row_index + grid_cols_padded, 0, dest_bytes_per_row);
			memset(Vy_max_grid + dest_row_index + grid_cols_padded, 0, dest_bytes_per_row);
			memset(Vy_grid + dest_row_index + grid_cols_padded, 0, dest_bytes_per_row);
			//memset(Fp_ywidth + dest_row_index + grid_cols_padded, 0, dest_bytes_per_row); PFU needs to be initialised before init subgrid

			if (Statesptr->routing == ON)
			{
				memset(route_V_ratio_per_sec_qy + dest_row_index + grid_cols_padded, 0, dest_bytes_per_row);
				memset(route_V_ratio_per_sec_qx + dest_row_index, 0, dest_bytes_per_row);
			}

			if (Statesptr->routing == ON)
			{
				InitRoutingStructure_row(route_V_ratio_per_sec_qx, route_V_ratio_per_sec_qy, j, not_last_row, grid_rows, grid_cols, grid_cols_padded, Arrptr);
			}
		}
	}
	printf("Init done\n");
	fflush(stdout);

#ifndef RESULT_CHECK

	memory_free_legacy(&Arrptr->Weir_Identx); //Weir_Identx converted to list, grid no longer needed
	memory_free_legacy(&Arrptr->Weir_Identy); //Weir_Identy converted to list, grid no longer needed

	// clean other memory that is no longer required after it has been copied/used
	memory_free_legacy(&Arrptr->FlowDir);
	memory_free_legacy(&Arrptr->RouteInt);
	memory_free_legacy(&Arrptr->H);
	//memory_free_legacy(& Arrptr->H);
	memory_free_legacy(&Arrptr->DEM);
	memory_free_legacy(&Arrptr->DEM);
	if (Arrptr->Manningsn != NULL)
		memory_free_legacy(&Arrptr->Manningsn);
	if (Arrptr->SGCManningsn != NULL)
		memory_free_legacy(&Arrptr->SGCManningsn);


	memory_free_legacy(&Arrptr->QxSGold);
	memory_free_legacy(&Arrptr->QySGold);

	memory_free_legacy(&Arrptr->SGCwidth);
	memory_free_legacy(&Arrptr->SGCz);
	memory_free_legacy(&Arrptr->SGCc);
	memory_free_legacy(&Arrptr->SGCbfH);
	memory_free_legacy(&Arrptr->SGCbfV);

	memory_free_legacy(&Arrptr->dx);
	memory_free_legacy(&Arrptr->dy);
	memory_free_legacy(&Arrptr->dA);
#endif
	Fast_IterateLoop(grid_cols, grid_rows, grid_cols_padded,
		h_grid, volume_grid,
		Qx_grid, Qy_grid, Qx_old_grid, Qy_old_grid,

		maxH_grid, maxHtm_grid, initHtm_grid, totalHtm_grid,
		maxVc_grid, maxVc_height_grid, maxHazard_grid,
		Vx_grid, Vy_grid, Vx_max_grid, Vy_max_grid,
		dem_grid,
		g_friction_sq_x_grid, g_friction_sq_y_grid,
		friction_x_grid, friction_y_grid,
		dx_col, dy_col, cell_area_col,
		Fp_xwidth,Fp_ywidth,

		sub_grid_layout_rows,
		sub_grid_state_rows,
		sub_grid_layout_blocks,
		sub_grid_state_blocks,

		SGC_BankFullHeight_grid,

		evap_time_series, &evap_grid, rain_time_series, rain_grid, dist_infil_grid, 

		wet_dry_bounds,
		ps_layout, boundaryCondition,
		weirs_weirs, weirs_bridges,
		route_dynamic_list,
		route_V_ratio_per_sec_qx, route_V_ratio_per_sec_qy,

		Super_linksptr,
		Fnameptr,
		Fptr,
		Locptr,
		Statesptr,
		Parptr,
		Solverptr,
		Damptr,
		SGCptr,
		tmp_thread_data,
#ifdef RESULT_CHECK
		Arrptr, BCptr,
		ChannelSegments,
		ChannelSegmentsVecPtr,
#endif

		verbose);

	if (verbose == ON) printf("Finished.\n\n");
}

// ITERATE THROUGH TIME STEPS
void Fast_MainStart(Fnames *Fnameptr, Files *Fptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr, Stage *Locptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, SGCprams *SGCptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, DamData *Damptr, int verbose)
{
	if (verbose == ON)
	{
		printf("\nStarting time steps: ");
		fflush(stdout);
	}
	Solverptr->itrn_time_now = Solverptr->itrn_time;

	// Populating Tstep variables prior to start of simulation
	if (Statesptr->adaptive_ts == ON)
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("adaptive mode\n\n");
		fflush(stdout);
	}
	else if (Statesptr->acceleration == ON)
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("acceleration mode\n\n");
		fflush(stdout);
	}
	else if (Statesptr->Roe == ON)
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("Roe mode\n\n");
		fflush(stdout);
	}
	else
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("non-adaptive mode\n\n");
		fflush(stdout);
	}
	if (Statesptr->SGC == ON)
	{
		// because the SGC model calculates the time step in UpdateH rather than during calcFPflow it needs to initalise 
		// SGCtmpTstep, which would usually be calculated in UpdateH
		Solverptr->Tstep = Solverptr->InitTstep;
		CalcT(Parptr, Solverptr, Arrptr);
		Solverptr->SGCtmpTstep = Solverptr->Tstep;
		if (verbose == ON) printf("SGC mode\n\n");
		fflush(stdout);
	}
	//	NUMERIC_TYPE init[16];
	//#pragma omp parallel for default(shared) schedule(static)
	//	for (int j = 0; j < 16; j++)
	//	{
	//		#pragma omp critical
	//		{
	//			init[j] = CBRT((NUMERIC_TYPE)j);
	//		}
	//
	//	}
	//	for (int j = 0; j < 16; j++)
	//	{
	//		printf("Init done: %" NUM_FMT"\n", init[j]);
	//	}

	Fast_MainInit(Fnameptr, Fptr, Statesptr, Parptr, Solverptr, BCptr, Locptr, ChannelSegments, Arrptr, SGCptr, ChannelSegmentsVecPtr, Damptr, verbose);
}
void InitSuperLinksStructure(const int grid_rows, const int grid_cols, const int grid_cols_padded, SuperGridLinksList * Super_linksptr, States *Statesptr, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr, Solver * Solverptr, Fnames * Fnameptr, int verbose)
{
	int n_links = 0, n_SGC_under_mask = 0;
	int i, j, p0, adjacent_to_chanmask;
	if (verbose == ON) printf("Creating SGC to 2D links list...  ");

	if (Statesptr->LinkListRead == ON) // read links from file
	{
		int i;
		NUMERIC_TYPE n;
		FILE *fp;
		fp = fopen_or_die(Fnameptr->LinkListfilename, "r", "Loading link list information", verbose);
		fscanf(fp, "%d", &n_links);
		fgetc(fp); // Retrieve closing EOL
		// Allocate memory
		AllocateSuperLinksMemory(n_links, Super_linksptr);
		// structure for file  2D_i, 2D_j, SGC_i, SGC_j, dx, w, n
		for (i = 0; i < n_links; i++)
		{
			fscanf(fp, "%i %i %i %i %" NUM_FMT" %" NUM_FMT" %" NUM_FMT"", &Super_linksptr->link_index_2D_i[i], &Super_linksptr->link_index_2D_j[i], &Super_linksptr->link_index_SGC_i[i], &Super_linksptr->link_index_SGC_j[i],
				&Super_linksptr->dx[i], &Super_linksptr->w[i], &n);
			Super_linksptr->gn2[i] = Solverptr->g*n*n;
			Super_linksptr->Qold[i] = C(0.0);
			Super_linksptr->link_index_2D[i] = Super_linksptr->link_index_2D_i[i] + Super_linksptr->link_index_2D_j[i] * grid_cols_padded;
			Super_linksptr->link_index_SGC[i] = Super_linksptr->link_index_SGC_i[i] + Super_linksptr->link_index_SGC_j[i] * grid_cols_padded;
			p0 = Super_linksptr->link_index_2D_i[i] + Super_linksptr->link_index_2D_j[i] *Parptr->xsz;
			Super_linksptr->DEM_z[i] = Arrptr->DEM[p0];
			p0 = Super_linksptr->link_index_SGC_i[i] + Super_linksptr->link_index_SGC_j[i] * Parptr->xsz;
			Super_linksptr->SGC_z[i] = Arrptr->SGCz[p0];
			Super_linksptr->SGC_bfH[i] = Arrptr->SGCbfH[p0];
		}
		fclose(fp);
		
	}
	else // calculate links in code using nearest neighbour
	{
		// Step -1 work out how many links and cells under mask will be needed 
		for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++) // loop through model domain
		{
			p0 = i + j * Parptr->xsz;       // cell index
			if (Arrptr->ChanMask[p0] > 0) // ChanMask cell
			{                             // Data Needed: location (x,y), bedz
				if (Arrptr->SGCwidth[p0] > C(0.0)) // identify if sub-grid channel cells and increment counter for memory allocation
				{
					n_SGC_under_mask = n_SGC_under_mask++;
				}
			}
			else if (Arrptr->SGCwidth[p0] == C(0.0))// 2D domain cell ... identify all locations on a ChanMask edge that do not include a sub-grid channel	 
			{    // Data Needed: location(x,y), z, n, w, dx
				adjacent_to_chanmask = 0;
				// Check north
				if (j > 0 && Arrptr->ChanMask[p0 - Parptr->xsz] > 0) adjacent_to_chanmask = 1;
				// Check south
				if (j + 1 < Parptr->ysz && Arrptr->ChanMask[p0 + Parptr->xsz] > 0) adjacent_to_chanmask = 1;
				// Check east
				if (i < Parptr->xsz - 1 && Arrptr->ChanMask[p0 + 1] > 0) adjacent_to_chanmask = 1;
				// Check west
				if (i > 0 && Arrptr->ChanMask[p0 - 1] > 0) adjacent_to_chanmask = 1;
				// Get the data needed
				if (adjacent_to_chanmask == 1)	n_links = n_links++;
			}
		}
		
		// Step 0: Create temp variables of the right size
		int* tmp_link_index_SGC_i = new int[n_SGC_under_mask]();
		int* tmp_link_index_SGC_j = new int[n_SGC_under_mask]();
		int* tmp_link_index_2D_i = new int[n_links]();
		int* tmp_link_index_2D_j = new int[n_links]();
		NUMERIC_TYPE* tmp_z_SGC = new NUMERIC_TYPE[n_SGC_under_mask]();
		NUMERIC_TYPE* tmp_z_2D = new NUMERIC_TYPE[n_links]();
		NUMERIC_TYPE* tmp_gn2 = new NUMERIC_TYPE[n_links]();
		NUMERIC_TYPE* tmp_w = new NUMERIC_TYPE[n_links]();
		NUMERIC_TYPE* tmp_dx = new NUMERIC_TYPE[n_links]();
		NUMERIC_TYPE* tmp_bfH = new NUMERIC_TYPE[n_SGC_under_mask]();

		// reset links counters to go again
		n_links = 0;
		n_SGC_under_mask = 0;

		// Step 1: Loop around ChanMask idetify locations of SGC and edges and gather data
		for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++) // loop through model domain
		{
			p0 = i + j * Parptr->xsz;       // cell index
			if (Arrptr->ChanMask[p0] > 0) // ChanMask cell
			{                             // Data Needed: location (x,y), bedz
				if (Arrptr->SGCwidth[p0] > C(0.0)) // identify if sub-grid channel cells // store location
				{
					tmp_link_index_SGC_i[n_SGC_under_mask] = i;
					tmp_link_index_SGC_j[n_SGC_under_mask] = j;
					tmp_z_SGC[n_SGC_under_mask] = Arrptr->SGCz[p0];
					tmp_bfH[n_SGC_under_mask] = Arrptr->SGCbfH[p0];
					n_SGC_under_mask = n_SGC_under_mask++;
				}
			}
			else if (Arrptr->SGCwidth[p0] == C(0.0))// 2D domain cell ... identify all locations on a ChanMask edge that do not include a sub-grid channel	 
			{    // Data Needed: location(x,y), z, n, w, dx
				adjacent_to_chanmask = 0;
				NUMERIC_TYPE dx_n = C(0.0), dx_e = C(0.0), dx_s = C(0.0), dx_w = C(0.0);
				NUMERIC_TYPE w_n = C(0.0), w_e = C(0.0), w_s = C(0.0), w_w = C(0.0);
				NUMERIC_TYPE max_dx, max_dy, max_wx, max_wy;
				// Check north
				if (j > 0 && Arrptr->ChanMask[p0 - Parptr->xsz] > 0)
				{
					adjacent_to_chanmask = 1; dx_n = Arrptr->dy[p0]; w_n = Arrptr->dx[p0];
				}
				// Check south
				if (j + 1 < Parptr->ysz && Arrptr->ChanMask[p0 + Parptr->xsz] > 0)
				{
					adjacent_to_chanmask = 1; dx_s = Arrptr->dy[p0]; w_s = Arrptr->dx[p0];
				}
				// Check east
				if (i < Parptr->xsz - 1 && Arrptr->ChanMask[p0 + 1] > 0)
				{
					adjacent_to_chanmask = 1; dx_e = Arrptr->dx[p0]; w_e = Arrptr->dy[p0];
				}
				// Check west
				if (i > 0 && Arrptr->ChanMask[p0 - 1] > 0)
				{
					adjacent_to_chanmask = 1; dx_w = Arrptr->dx[p0]; w_w = Arrptr->dy[p0];
				}
				// Get the data needed
				if (adjacent_to_chanmask == 1)
				{
					// store location
					tmp_link_index_2D_i[n_links] = i;
					tmp_link_index_2D_j[n_links] = j;
					// get variables
					tmp_z_2D[n_links] = Arrptr->DEM[p0]; // DEM z
					if (Arrptr->Manningsn != NULL)   tmp_gn2[n_links] = (Arrptr->Manningsn[p0] * Arrptr->Manningsn[p0]); // Manning's n
					else tmp_gn2[n_links] = Parptr->FPn * Parptr->FPn;
					tmp_gn2[n_links] = tmp_gn2[n_links] * Solverptr->g;
					tmp_w[n_links] = FMAX(w_n, w_s) + FMAX(w_e, w_w);  // flow width (max half of cell edge)
					// calculate dx to use
					if (FMAX(dx_n, dx_s) > C(0.0) && FMAX(dx_e, dx_w) > C(0.0))	tmp_dx[n_links] = sqrt(FMAX(dx_n, dx_s) * FMAX(dx_n, dx_s) + FMAX(dx_e, dx_w) * FMAX(dx_e, dx_w)); // fow from x and y
					else if (FMAX(dx_n, dx_s) > C(0.0)) tmp_dx[n_links] = FMAX(dx_n, dx_s);
					else if (FMAX(dx_e, dx_w) > C(0.0)) tmp_dx[n_links] = FMAX(dx_e, dx_w);
					else printf("Error in SGC to DEM link function");
					n_links = n_links++;
				}
			}
		}

		// Step 2: Allocate memory
		AllocateSuperLinksMemory(n_links, Super_linksptr);

		// Step 3: for each edge identify closest sug_grid cell for linking list
#pragma omp for private (j)
		for (i = 0; i < n_links; i++)
		{
			// Copy 2D info to structure
			Super_linksptr->link_index_2D_i[i] = tmp_link_index_2D_i[i];
			Super_linksptr->link_index_2D_j[i] = tmp_link_index_2D_j[i];
			Super_linksptr->link_index_2D[i] = tmp_link_index_2D_i[i] + tmp_link_index_2D_j[i] * grid_cols_padded;
			Super_linksptr->DEM_z[i] = tmp_z_2D[i];
			Super_linksptr->dx[i] = tmp_dx[i];
			Super_linksptr->w[i] = tmp_w[i];
			Super_linksptr->gn2[i] = tmp_gn2[i];
			Super_linksptr->Qold[i] = C(0.0);

			// now work out closest sub-grid cell
			int closest_i, closest_j, closest_line;
			//NUMERIC_TYPE dist, dist2;
			int dist, dist2, tmp1, tmp2;
			//dist = sqrt(pow((NUMERIC_TYPE)tmp_link_index_2D_i[i] - (NUMERIC_TYPE)tmp_link_index_SGC_i[0], 2.0) + pow((NUMERIC_TYPE)tmp_link_index_2D_j[i] - (NUMERIC_TYPE)tmp_link_index_SGC_j[0], 2.0));
			dist = 2000000000; // set to big number
			closest_i = tmp_link_index_SGC_i[0]; closest_j = tmp_link_index_SGC_j[0], closest_line = 0;
			for (j = 0; j < n_SGC_under_mask; j++) // loop through sub-grid cells
			{
				//dist2 = sqrt(pow((NUMERIC_TYPE)tmp_link_index_2D_i[i] - (NUMERIC_TYPE)tmp_link_index_SGC_i[j], 2.0) + pow((NUMERIC_TYPE)tmp_link_index_2D_j[i] - (NUMERIC_TYPE)tmp_link_index_SGC_j[j], 2.0));
				tmp1 = tmp_link_index_2D_i[i] - tmp_link_index_SGC_i[j];
				tmp2 = tmp_link_index_2D_j[i] - tmp_link_index_SGC_j[j];
				tmp1 *= tmp1; tmp2 *= tmp2;
				dist2 = tmp1 + tmp2;
				if (dist > dist2)
				{
					closest_line = j;
					dist = dist2;
				}
			}
			closest_i = tmp_link_index_SGC_i[closest_line];
			closest_j = tmp_link_index_SGC_j[closest_line];

			Super_linksptr->link_index_SGC_i[i] = closest_i;
			Super_linksptr->link_index_SGC_j[i] = closest_j;
			Super_linksptr->link_index_SGC[i] = closest_i + closest_j * grid_cols_padded;
			Super_linksptr->SGC_z[i] = tmp_z_SGC[closest_line];
			Super_linksptr->SGC_bfH[i] = tmp_bfH[closest_line];
		}
		if (verbose == ON)
		{
			printf("%i links created.\n", n_links);
			printf("Connecting to %i sub-grid channels. Done.\n", n_SGC_under_mask);
		}
	}
}
void AllocateSuperLinksMemory(int n_links, SuperGridLinksList * Super_linksptr)
{
	// Step 2: Allocate memory
	Super_linksptr->num_links = n_links;
	Super_linksptr->DEM_z = (NUMERIC_TYPE*)memory_allocate(n_links * sizeof(NUMERIC_TYPE));
	Super_linksptr->SGC_z = (NUMERIC_TYPE*)memory_allocate(n_links * sizeof(NUMERIC_TYPE));
	Super_linksptr->SGC_bfH = (NUMERIC_TYPE*)memory_allocate(n_links * sizeof(NUMERIC_TYPE));
	Super_linksptr->dx = (NUMERIC_TYPE*)memory_allocate(n_links * sizeof(NUMERIC_TYPE));
	Super_linksptr->w = (NUMERIC_TYPE*)memory_allocate(n_links * sizeof(NUMERIC_TYPE));
	Super_linksptr->gn2 = (NUMERIC_TYPE*)memory_allocate(n_links * sizeof(NUMERIC_TYPE));
	Super_linksptr->Qold = (NUMERIC_TYPE*)memory_allocate(n_links * sizeof(NUMERIC_TYPE));
	Super_linksptr->link_index_SGC_i = (int*)memory_allocate(n_links * sizeof(int));
	Super_linksptr->link_index_SGC_j = (int*)memory_allocate(n_links * sizeof(int));
	Super_linksptr->link_index_2D_i = (int*)memory_allocate(n_links * sizeof(int));
	Super_linksptr->link_index_2D_j = (int*)memory_allocate(n_links * sizeof(int));
	Super_linksptr->link_index_2D = (int*)memory_allocate(n_links * sizeof(int));
	Super_linksptr->link_index_SGC = (int*)memory_allocate(n_links * sizeof(int));
}
