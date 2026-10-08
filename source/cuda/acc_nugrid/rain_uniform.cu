//#include "rain.cuh"
//#include "../../rain/rain.h"
//#include "unifiedallocator.cuh"
//#include "BLOCK_VAR_MACROS.cuh"

#include "rain_uniform.cuh"

__global__ 
void lis::cuda::acc_nugrid::update_uniform_rain
(
    AssembledSolution    d_assem_sol,
    Pars pars,
    NUMERIC_TYPE                 dt,
    NUMERIC_TYPE  rain_rate
)
{
	index_1D t_idx = threadIdx.x;
	index_1D idx = blockIdx.x * blockDim.x + t_idx;

    if (idx >= d_assem_sol.length) return;

    NUMERIC_TYPE cell_rain;
    
    cell_rain = rain_rate * dt; //rate for depth, not area.  Has to be inside loop due to negative rain legacy support.

    NUMERIC_TYPE Z = d_assem_sol.z0[idx];
    if (FABS(Z - pars.nodata_elevation) < C(1e-6)) {
        //continue; do nothing
    }
    else {

        //check for -ve depths (legacy support for pre-evap code when rainfall could be negative)
        //if (h0 < 0)
        //{
        //    cell_rain += h0;
        //    h0 = 0;
        //}
        
        d_assem_sol.h[idx] += cell_rain;

//        loc_rainfall_total += cell_rain * Parptr->dA; // mass balance for local cell (cumulative)	
    }

}