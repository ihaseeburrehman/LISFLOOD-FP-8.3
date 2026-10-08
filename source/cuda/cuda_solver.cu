#include "cuda_solver.cuh"
#include "cuda_atomic.cuh"

__constant__ Geometry lis::cuda::geometry;
__constant__ Geometry lis::cuda::rain_geometry;
__constant__ int lis::cuda::pitch;
__constant__ lis::PhysicalParams lis::cuda::physical_params;
__constant__ lis::SolverParams lis::cuda::solver_params;
__constant__ lis::cuda::BoundaryConditions lis::cuda::boundaries;
__constant__ lis::SamplePoints lis::cuda::sample_points;
__managed__ NUMERIC_TYPE lis::cuda::dt;
__managed__ int lis::cuda::sample_buf_idx = 0;

const dim3 lis::cuda::block_size(CUDA_BLOCK_SIZE_X, CUDA_BLOCK_SIZE_Y);

__device__ void lis::cuda::update_mass_stats_x
(
	MassStats* stats,
	NUMERIC_TYPE FH,
	int i,
	int j,
	NUMERIC_TYPE time_stage_fraction
)
{
	if (i == 0)
	{
		if (FH > C(0.0))
		{
			atomicAdd(&(stats->in), time_stage_fraction*FH*cuda::geometry.dy);
		}
		else
		{
			atomicAdd(&(stats->out), -time_stage_fraction*FH*cuda::geometry.dy);
		}
	}
	else if (i == cuda::geometry.xsz)
	{
		if (FH < C(0.0))
		{
			atomicAdd(&(stats->in), -time_stage_fraction*FH*cuda::geometry.dy);
		}
		else
		{
			atomicAdd(&(stats->out), time_stage_fraction*FH*cuda::geometry.dy);
		}
	}
}

__device__ void lis::cuda::update_mass_stats_y
(
	MassStats* stats,
	NUMERIC_TYPE FH,
	int i,
	int j,
	NUMERIC_TYPE time_stage_fraction
)
{
	if (j == 0)
	{
		if (FH < C(0.0))
		{
			atomicAdd(&(stats->in), -time_stage_fraction*FH*cuda::geometry.dx);
		}
		else
		{
			atomicAdd(&(stats->out), time_stage_fraction*FH*cuda::geometry.dx);
		}
	}
	else if (j == cuda::geometry.ysz)
	{
		if (FH > C(0.0))
		{
			atomicAdd(&(stats->in), time_stage_fraction*FH*cuda::geometry.dx);
		}
		else
		{
			atomicAdd(&(stats->out), -time_stage_fraction*FH*cuda::geometry.dx);
		}
	}
}
