#include "rain.h"
#include "../lisflood2/file_tool.h"
#include "../utility.h"

static const NUMERIC_TYPE epsilon = C(1e-12);

template<class Allocator>
DynamicRain<Allocator>::DynamicRain
(
    const char* filename,
    int verbose,
    const Allocator& allocator
)
:
allocator_(allocator)
{
	enabled_ = (strlen(filename) > 0);
    if (!enabled_) return;

    read_file_netCDF_start(filename, "rainfall_depth", &netcdf_);
    netcdf_.data = allocator_.allocate(netcdf_.xlen * netcdf_.ylen);

    geometry_.xsz = netcdf_.xlen;
    geometry_.ysz = netcdf_.ylen;
    geometry_.dx = netcdf_.xs[1] - netcdf_.xs[0];
    geometry_.dy = netcdf_.ys[0] - netcdf_.ys[1];
    geometry_.blx = netcdf_.xs[0] - geometry_.dx/C(2.0);
    geometry_.tly = netcdf_.ys[0] + geometry_.dy/C(2.0);
    geometry_.bly = geometry_.tly - geometry_.ysz*geometry_.dy;

    update_time(C(0.0));
}

template<class Allocator>
bool DynamicRain<Allocator>::has_same_origin
(
    Pars* Parptr
)
{
    return FABS(Parptr->blx - geometry_.blx) < epsilon &&
        FABS(Parptr->bly - geometry_.bly) < epsilon;
}

template<class Allocator>
bool DynamicRain<Allocator>::is_tile_size_multiple_of_grid
(
    Pars* Parptr
)
{
    NUMERIC_TYPE dx_ratio = geometry_.dx / Parptr->dx;
    NUMERIC_TYPE dy_ratio = geometry_.dy / Parptr->dy;

    return FABS(round(dx_ratio) - dx_ratio) <= epsilon &&
            FABS(round(dy_ratio) - dy_ratio) <= epsilon;
}

template<class Allocator>
void DynamicRain<Allocator>::update_time
(
    NUMERIC_TYPE t
)
{
    if (!enabled_) return;

    if (read_file_netCDF(&netcdf_, t / C(3600.0) /* s to hr */))
    {
        netcdf_.dt *= C(3600.0); // hr to s

        // convert from rainfall accumulated over dt (mm) to rainfall rate (m/s)
        for (int i=0; i<geometry_.xsz*geometry_.ysz; i++)
        {
            netcdf_.data[i] /= netcdf_.dt;
            netcdf_.data[i] /= C(1000.0); // mm to m
        }
    }
}

template<class Allocator>
bool DynamicRain<Allocator>::update_time_SGC
(
    NUMERIC_TYPE t
)
{
    if (!enabled_) return (false);

	bool time_updated;
	time_updated = read_file_netCDF(&netcdf_, t / C(3600.0) /* s to hr */);
    if (time_updated)
    {
        netcdf_.dt *= C(3600.0); // hr to s

        // convert from rainfall accumulated over dt (mm) to rainfall rate (m/s)
        for (int i=0; i<geometry_.xsz*geometry_.ysz; i++)
        {
            netcdf_.data[i] /= netcdf_.dt;
            netcdf_.data[i] /= C(1000.0); // mm to m
        }
    }
    return time_updated;
}


template<class Allocator>
NUMERIC_TYPE DynamicRain<Allocator>::rate_at_cell
(
    Pars* Parptr,
    int i,
    int j
)
{
    if (!enabled_) return C(0.0);

    int tile_i = i * Parptr->dx / geometry_.dx;

    NUMERIC_TYPE top_gap = Parptr->tly - geometry_.tly;
    int tile_j = (j - top_gap/Parptr->dy) * Parptr->dy / geometry_.dy;

    if (tile_i < geometry_.xsz && tile_j >= 0)
    {
        return netcdf_.data[tile_j*geometry_.xsz + tile_i];
    }
    else
    {
        return C(0.0);
    }
}


template<class Allocator>
NUMERIC_TYPE DynamicRain<Allocator>::rate_at_cell_SGC
(
    int i,
    int j,
    const NUMERIC_TYPE dx_col, 
    const NUMERIC_TYPE dy_col,
    NUMERIC_TYPE tly
)
{
    if (!enabled_) return C(0.0);

    int tile_i = i * dx_col / geometry_.dx;

    NUMERIC_TYPE top_gap = tly - geometry_.tly;
    int tile_j = (j - top_gap/dy_col) * dy_col / geometry_.dy;

    if (tile_i < geometry_.xsz && tile_j >= 0)
    {
        return netcdf_.data[tile_j*geometry_.xsz + tile_i];
    }
    else
    {
        return C(0.0);
    }
}

template<class Allocator>
void DynamicRain<Allocator>::update_H
(
    Pars *Parptr,
    Solver *Solverptr,
    Arrays *Arrptr
)
{
#if _NETCDF == 1
    if (!enabled_) return;

    update_time(Solverptr->t);
    NUMERIC_TYPE total_rain_mass = C(0.0);
#pragma omp parallel for reduction (+:total_rain_mass)
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
            if (FABS(Z - Parptr->nodata_elevation) < 1e-6) continue;

			NUMERIC_TYPE& H = Arrptr->H[j*Parptr->xsz + i];
            NUMERIC_TYPE inc = rate_at_cell(Parptr, i, j) * Solverptr->Tstep;
            H += inc;
            total_rain_mass += inc*Parptr->dA;
        }
    }
#endif
    Parptr->RainTotalLoss += total_rain_mass;
}

template<class Allocator>
void DynamicRain<Allocator>::update_rain_grid_SGM
(
    const NUMERIC_TYPE curr_time,
    NUMERIC_TYPE * rain_grid,
    const NUMERIC_TYPE * dem_grid,
    WetDryRowBound * wet_dry_bounds,
    const NUMERIC_TYPE *dx_col, 
    const NUMERIC_TYPE *dy_col,
    const NUMERIC_TYPE tly,
    const int grid_rows,
    const int grid_cols_padded,
    const NUMERIC_TYPE *cell_area_col
)
{
#if _NETCDF == 1
  if (!enabled_) return;

  //update_time(curr_time); 
  if (update_time_SGC(curr_time)); // if the file had to be read update rain_grid
  {	
	for (int j = 0; j < grid_rows; j++)
	{
		// update wet_dry_bounds as all dem cells will now be wet
		wet_dry_bounds->fp_vol[j] = wet_dry_bounds->dem_data[j];
		const int row_start = wet_dry_bounds->dem_data[j].start;
		const int row_end = wet_dry_bounds->dem_data[j].end;
		int grid_row_index = j * grid_cols_padded;
	
#ifdef __INTEL_COMPILER
	__assume_aligned(rain_grid, 64);
#endif
#pragma ivdep
		for (int i = row_start; i < row_end; i++)
		{
			int grid_cell_index = grid_row_index + i;
			const NUMERIC_TYPE row_dx = dx_col[j];
			const NUMERIC_TYPE row_dy = dy_col[j];
			if (dem_grid[grid_cell_index] != DEM_NO_DATA)
			{
				rain_grid[grid_cell_index] = rate_at_cell_SGC(i, j,row_dx, row_dy, tly)*cell_area_col[j];	
			}
		}
	}
  }
#endif
  }


template<class Allocator>
DynamicRain<Allocator>::~DynamicRain()
{
    if (!enabled_) return;
    free(netcdf_.xs);
    free(netcdf_.ys);
    free(netcdf_.times);
    allocator_.deallocate(netcdf_.data, netcdf_.xlen * netcdf_.ylen);
    CloseNetCDF(netcdf_.ncid);
}
