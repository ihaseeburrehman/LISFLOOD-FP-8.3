#include "cuda_rain.cuh"
#include "../rain/rain.h"
#include "cuda_unifiedallocator.cuh"

#if _NETCDF == 1
	#include "../rain/rain.tpp"
#else
	#include "../rain/rain_stub.tpp"
#endif

template class DynamicRain<lis::cuda::UnifiedAllocator<NUMERIC_TYPE>>;

__device__ 
NUMERIC_TYPE lis::cuda::rate_at_cell
(
	NUMERIC_TYPE* rain_data,
	int i,
	int j
)
{
    int tile_i = (i-1) * cuda::geometry.dx / cuda::rain_geometry.dx;

    NUMERIC_TYPE top_gap = cuda::geometry.tly - cuda::rain_geometry.tly;
    int tile_j = (j-1 - top_gap/cuda::geometry.dy) *
		cuda::geometry.dy / cuda::rain_geometry.dy;

    if (tile_i < cuda::rain_geometry.xsz && tile_j >= 0)
    {
        return rain_data[tile_j*cuda::rain_geometry.xsz + tile_i];
    }
    else
    {
        return C(0.0);
    }
}

__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
lis::cuda::DynamicRain::update
(
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* rain_data
)
{
#if _NETCDF == 1
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j+1; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i+1; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z = DEM[j*cuda::pitch + i];
            if (FABS(Z - cuda::solver_params.nodata_elevation) < C(1e-6)) continue;

			NUMERIC_TYPE& Hval = H[j*cuda::pitch + i];
            Hval += rate_at_cell(rain_data, i, j) * cuda::dt;
		}
	}
#endif
}


__device__ 
NUMERIC_TYPE lis::cuda::rate_at_cell_ACC
(
	NUMERIC_TYPE* rain_data,
	int i,
	int j
)
{
    int tile_i = (i) * cuda::geometry.dx / cuda::rain_geometry.dx;

    NUMERIC_TYPE top_gap = cuda::geometry.tly - cuda::rain_geometry.tly;
    int tile_j = (j - top_gap/cuda::geometry.dy) *
		cuda::geometry.dy / cuda::rain_geometry.dy;

    if (tile_i < cuda::rain_geometry.xsz && tile_j >= 0)
    {
        return rain_data[tile_j*cuda::rain_geometry.xsz + tile_i];
    }
    else
    {
        return C(0.0);
    }
}

__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
lis::cuda::DynamicRain::updateACC
(
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* rain_data
)
{
#if _NETCDF == 1
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j; j<cuda::geometry.ysz; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i; i<cuda::geometry.xsz; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z = DEM[j*cuda::geometry.xsz + i];
            if (FABS(Z - cuda::solver_params.nodata_elevation) < C(1e-6)) continue;

			NUMERIC_TYPE& Hval = H[j*cuda::geometry.xsz + i];
            Hval += rate_at_cell_ACC(rain_data, i, j) * cuda::dt;
		}
	}
#endif
}