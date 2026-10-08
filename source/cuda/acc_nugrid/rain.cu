#include "rain.cuh"
#include "../../rain/rain.h"
#include "unifiedallocator.cuh"
#include "BLOCK_VAR_MACROS.cuh"

#if _NETCDF == 1
	#include "../../rain/rain.tpp"
#else
	#include "../../rain/rain_stub.tpp"
#endif

template class DynamicRain<UnifiedAllocator<NUMERIC_TYPE>>;

__device__ 
NUMERIC_TYPE lis::cuda::acc_nugrid::rate_at_cell
(
	NUMERIC_TYPE* rain_data,
	int i,
	int j,
    NUMERIC_TYPE dx,
    Pars pars,
    int xsz_rain,
    NUMERIC_TYPE tly_rain,
    NUMERIC_TYPE dx_rain,
    NUMERIC_TYPE dy_rain
)
{
    int tile_i = (i) * dx / dx_rain;

    NUMERIC_TYPE top_gap = pars.tly - tly_rain;
    int tile_j = (j - top_gap/ dx) * dx / dy_rain;

    if (tile_i < xsz_rain && tile_j >= 0)
    {
        return rain_data[tile_j*xsz_rain + tile_i];
    }
    else
    {
        return C(0.0);
    }
}

__global__ 
void lis::cuda::acc_nugrid::update_rain
(
    AssembledSolution    d_assem_sol,
	NUMERIC_TYPE* rain_data,
    int xsz_rain,
    NUMERIC_TYPE tly_rain,
    NUMERIC_TYPE dx_rain,
    NUMERIC_TYPE dy_rain,
    Pars pars,
    Solver solver
)
{
#if _NETCDF == 1
	index_1D t_idx = threadIdx.x;
	index_1D idx = blockIdx.x * blockDim.x + t_idx;

    if (idx >= d_assem_sol.length) return;

    NUMERIC_TYPE dx = pars.dx * (1 << (solver.L - d_assem_sol.levels[idx]));

    int level = d_assem_sol.levels[idx];

    MortonCode code = d_assem_sol.act_idcs[idx] - get_lvl_idx(level); // MKS: code in level

    int i = compact(code);
    int j = compact(code >> 1);

    NUMERIC_TYPE Z = d_assem_sol.z0[idx];
    if (FABS(Z - pars.nodata_elevation) < C(1e-6)) {
        //continue; do nothing
    }
    else {
        d_assem_sol.h[idx] += rate_at_cell(rain_data, i, j, dx, pars, xsz_rain, tly_rain, dx_rain, dy_rain) * solver.Tstep;
    }
#endif
}

__global__
void lis::cuda::acc_nugrid::drain_nodata_water
(
    AssembledSolution    d_assem_sol,
    Pars pars,
    Solver solver
)
{

    index_1D t_idx = threadIdx.x;
    index_1D idx = blockIdx.x * blockDim.x + t_idx;

    if (idx >= d_assem_sol.length) return;

    NUMERIC_TYPE dx = pars.dx * (1 << (solver.L - d_assem_sol.levels[idx]));

    int level = d_assem_sol.levels[idx];

    MortonCode code = d_assem_sol.act_idcs[idx] - get_lvl_idx(level); // MKS: code in level

    int i = compact(code);
    int j = compact(code >> 1);

    NUMERIC_TYPE Z = d_assem_sol.z0[idx];
    if (FABS(Z - pars.nodata_elevation) < C(1e-6)) {

        d_assem_sol.h[idx] = C(0.0);
    }

}

