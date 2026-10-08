#include "DataTypes.h"

void AllocateWetDryRowBound(int row_count, int block_count, WetDryRowBound * wet_dry_bounds)
{
	wet_dry_bounds->fp_h = (IndexRange*)memory_allocate(row_count * sizeof(IndexRange));
	wet_dry_bounds->fp_h_prev = (IndexRange*)memory_allocate(row_count * sizeof(IndexRange));
	wet_dry_bounds->fp_vol = (IndexRange*)memory_allocate(row_count * sizeof(IndexRange));
	wet_dry_bounds->dem_data = (IndexRange*)memory_allocate(row_count * sizeof(IndexRange));


	wet_dry_bounds->block_count = block_count;
	wet_dry_bounds->block_row_bounds = (IndexRange*)memory_allocate(sizeof(IndexRange*) * block_count);

	for (int block_index = 0; block_index < block_count; block_index++)
	{
		wet_dry_bounds->block_row_bounds[block_index].start = -1;
		wet_dry_bounds->block_row_bounds[block_index].end = -1;
	}
}

void AllocateSubGridCellInfo(int cell_count, SubGridCellInfo * sub_grid_cell_info)
{
	sub_grid_cell_info->cell_count = cell_count;

	sub_grid_cell_info->sg_cell_x = (int*)memory_allocate(cell_count * sizeof(int));
	sub_grid_cell_info->sg_cell_y = (int*)memory_allocate(cell_count * sizeof(int));
	sub_grid_cell_info->sg_cell_grid_index_lookup = (int*)memory_allocate(cell_count * sizeof(int));

	sub_grid_cell_info->sg_cell_cell_area = (NUMERIC_TYPE*)memory_allocate(cell_count * sizeof(NUMERIC_TYPE));
	sub_grid_cell_info->sg_cell_dem = (NUMERIC_TYPE*)memory_allocate(cell_count * sizeof(NUMERIC_TYPE));

	sub_grid_cell_info->sg_cell_SGC_width = (NUMERIC_TYPE*)memory_allocate(cell_count * sizeof(NUMERIC_TYPE));
	sub_grid_cell_info->sg_cell_SGC_c = (NUMERIC_TYPE*)memory_allocate(cell_count * sizeof(NUMERIC_TYPE));
	sub_grid_cell_info->sg_cell_SGC_BankFullHeight = (NUMERIC_TYPE*)memory_allocate(cell_count * sizeof(NUMERIC_TYPE));
	sub_grid_cell_info->sg_cell_SGC_BankFullVolume = (NUMERIC_TYPE*)memory_allocate(cell_count * sizeof(NUMERIC_TYPE));

	sub_grid_cell_info->sg_cell_SGC_group = (int*)memory_allocate(cell_count * sizeof(int));
	sub_grid_cell_info->sg_cell_SGC_is_large = (int*)memory_allocate(cell_count * sizeof(int));
}

void ZeroSubGridCellInfo(SubGridCellInfo * sub_grid_cell_info, int cell_index)
{
	sub_grid_cell_info->sg_cell_x[cell_index] = -1;
	sub_grid_cell_info->sg_cell_y[cell_index] = -1;
	sub_grid_cell_info->sg_cell_grid_index_lookup[cell_index] = -1;

	sub_grid_cell_info->sg_cell_cell_area[cell_index] = C(-1.0);
	sub_grid_cell_info->sg_cell_dem[cell_index] = C(-1.0);

	sub_grid_cell_info->sg_cell_SGC_width[cell_index] = C(-1.0);
	sub_grid_cell_info->sg_cell_SGC_BankFullHeight[cell_index] = C(-1.0);
	sub_grid_cell_info->sg_cell_SGC_BankFullVolume[cell_index] = C(-1.0);
	sub_grid_cell_info->sg_cell_SGC_c[cell_index] = C(-1.0);

	sub_grid_cell_info->sg_cell_SGC_group[cell_index] = -1;
	sub_grid_cell_info->sg_cell_SGC_is_large[cell_index] = -1;
}

void AllocateWaterSource(int count, WaterSource * waterSource)
{
	waterSource->count = count;

	waterSource->Ident = (ESourceType*)memory_allocate(count * sizeof(ESourceType));
	waterSource->Val = (NUMERIC_TYPE*)memory_allocate(count * sizeof(NUMERIC_TYPE));
	waterSource->timeSeries = (TimeSeries**)memory_allocate(count * sizeof(TimeSeries*));
	waterSource->Q_FP_old = (NUMERIC_TYPE*)memory_allocate(count * sizeof(NUMERIC_TYPE));
	waterSource->Q_SG_old = (NUMERIC_TYPE*)memory_allocate(count * sizeof(NUMERIC_TYPE));
	waterSource->g_friction_squared_FP = (NUMERIC_TYPE*)memory_allocate(count * sizeof(NUMERIC_TYPE));
	waterSource->g_friction_squared_SG = (NUMERIC_TYPE*)memory_allocate(count * sizeof(NUMERIC_TYPE));

	AllocateSubGridCellInfo(count, &waterSource->ws_cell);
}


void AllocateWeir(int count, WeirLayout * weirs)
{
	weirs->weir_index_qx = (int*)memory_allocate(count * sizeof(int));
	weirs->weir_index_qy = (int*)memory_allocate(count * sizeof(int));

	int weir_count = weirs->weir_count;

	weirs->Weir_Q_old_SG = (NUMERIC_TYPE*)memory_allocate(weir_count * sizeof(NUMERIC_TYPE));
	weirs->Weir_grid_index = (int*)memory_allocate(weir_count * sizeof(int));
	weirs->Weir_g_friction_sq = (NUMERIC_TYPE*)memory_allocate(weir_count * sizeof(NUMERIC_TYPE));

	weirs->Weir_pair_stream_flow_index = (int*)memory_allocate(2 * weir_count * sizeof(int));
	AllocateSubGridCellInfo(2 * weir_count, &weirs->cell_pair);
}


void AllocateRoutingDynamicList(int rows, int grid_cols_padded, RouteDynamicList * route_dynamic_list)
{
	route_dynamic_list->row_route_qx_count = (int*)memory_allocate(rows * sizeof(int));
	route_dynamic_list->row_route_qy_count = (int*)memory_allocate(rows * sizeof(int));
	route_dynamic_list->route_list_i_lookup_qx = (int*)memory_allocate(grid_cols_padded * rows * sizeof(int));
	route_dynamic_list->route_list_i_lookup_qy = (int*)memory_allocate(grid_cols_padded * rows * sizeof(int));
}