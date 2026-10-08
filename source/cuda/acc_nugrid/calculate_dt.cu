#include "calculate_dt.cuh"

__global__ void lis::cuda::acc_nugrid::calculate_dt
(
    AssembledSolution    d_assem_sol,
    Pars pars,
    Solver solver,
    NUMERIC_TYPE*                d_dt_CFL
)
{
    
    //int idx = blockIdx.x * blockDim.x + threadIdx.x;

    index_1D t_idx = threadIdx.x;
    index_1D idx = blockIdx.x * blockDim.x + t_idx;

//    if (idx >= d_assem_sol.length) return;

    if (idx < d_assem_sol.length)
    {
        
        d_dt_CFL[idx] = solver.InitTstep;
    }

    //d_dt_CFL[idx] = InitTstep;

    NUMERIC_TYPE dt_loc;

    NUMERIC_TYPE h = d_assem_sol.h[idx];

    NUMERIC_TYPE dx = pars.dx * (1 << (solver.L - d_assem_sol.levels[idx]));

    bool below_depth = (h < solver.DepthThresh);

    if (below_depth)
    {
        //return;
        d_dt_CFL[idx] = solver.InitTstep;
    }
    else
    {
        dt_loc = solver.cfl * dx / SQRT(solver.g * h);

        d_dt_CFL[idx] = FMIN(solver.Tstep, dt_loc);
    }


}
