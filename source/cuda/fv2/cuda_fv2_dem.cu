#include "cuda_fv2_dem.cuh"
#include "cuda_dem.cuh"
#include "cuda_geometry.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"
#include "io.h"

namespace lis
{
namespace cuda
{
namespace fv2
{

__global__ void initialise_Zstar_x
(
	NUMERIC_TYPE* Zstar_x,
	DeviceTopography DEM
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j+1; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z_neg = DEM._0[j*cuda::pitch + i];
			NUMERIC_TYPE Z_pos = DEM._0[j*cuda::pitch + i+1];

			Zstar_x[j*cuda::pitch + i] = FMAX(Z_neg, Z_pos);
		}
	}
}

__global__ void initialise_Zstar_y
(
	NUMERIC_TYPE* Zstar_y,
	DeviceTopography DEM
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i+1; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z_neg = DEM._0[(j+1)*cuda::pitch + i];
			NUMERIC_TYPE Z_pos = DEM._0[j*cuda::pitch + i];

			Zstar_y[j*cuda::pitch + i] = FMAX(Z_neg, Z_pos);
		}
	}
}

}
}
}

lis::cuda::fv2::Topography::Topography
(
	const char* filename,
	Geometry& geometry,
	int& pitch,
	int& offset,
	NUMERIC_TYPE nodata_elevation,
	int acceleration, //// added
	int verbose
)
{
	_0 = cuda::Topography::load(filename, geometry, pitch, offset,
			nodata_elevation, acceleration, verbose); //// modified

	clamp_boundary_values(geometry, pitch);
}




void lis::cuda::fv2::Topography::clamp_boundary_values
(
	Geometry& geometry,
	int pitch
)
{
	for (int j=1; j<geometry.ysz+1; j++)
	{
		// west
		{
			const int i=0;
			_0[j * pitch + i] = _0[j * pitch + i + 1];	
		}

		// east
		{
			const int i=geometry.xsz+1;
			_0[j * pitch + i] = _0[j * pitch + i - 1];
		}
	}

	for (int i=1; i<geometry.xsz+1; i++)
	{
		// north
		{
			const int j=0;
			_0[j * pitch + i] = _0[(j + 1) * pitch + i];
		}

		// south
		{
			const int j=geometry.ysz+1;
			_0[j * pitch + i] = _0[(j - 1) * pitch + i];
		}
	}
}

lis::cuda::fv2::Topography::~Topography()
{
	delete[] _0;
}

void lis::cuda::fv2::DeviceTopography::initialise
(
	Topography& DEM,
	Geometry& geometry
)
{
	_0 = cuda::GhostRaster::allocate_device(geometry);
	cuda::GhostRaster::copy(_0, DEM._0, geometry);

	Zstar_x = cuda::GhostRaster::allocate_device(geometry);
	Zstar_y = cuda::GhostRaster::allocate_device(geometry);
	initialise_Zstar_x();
	initialise_Zstar_y();
}


__device__
NUMERIC_TYPE lis::cuda::fv2::DeviceTopography::Zdagger
(
	NUMERIC_TYPE ETA,
	NUMERIC_TYPE Zstar
)
{
	return Zstar - FMAX(C(0.0), -(ETA - Zstar));
}

void lis::cuda::fv2::DeviceTopography::initialise_Zstar_x()
{
	lis::cuda::fv2::initialise_Zstar_x<<<1, cuda::block_size>>>(Zstar_x, *this);
}

void lis::cuda::fv2::DeviceTopography::initialise_Zstar_y()
{
	lis::cuda::fv2::initialise_Zstar_y<<<1, cuda::block_size>>>(Zstar_y, *this);
}

void lis::cuda::fv2::DeviceTopography::free_device()
{
	cuda::free_device(_0);
}
