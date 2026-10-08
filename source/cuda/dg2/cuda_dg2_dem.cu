#include "cuda_dg2_dem.cuh"
#include "cuda_dem.cuh"
#include "cuda_geometry.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"
#include "io.h"

namespace lis
{
namespace cuda
{
namespace dg2
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
			NUMERIC_TYPE Z_neg = DEM.neg_x(j*cuda::pitch + i);

			NUMERIC_TYPE Z_pos = DEM.pos_x(j*cuda::pitch + i+1);

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
			NUMERIC_TYPE Z_neg = DEM.neg_y((j+1)*cuda::pitch + i);

			NUMERIC_TYPE Z_pos = DEM.pos_y(j*cuda::pitch + i);

			Zstar_y[j*cuda::pitch + i] = FMAX(Z_neg, Z_pos);
		}
	}
}

}
}
}

lis::cuda::dg2::Topography::Topography
(
	const char* filename,
	Geometry& geometry,
	int& pitch,
	int& offset,
	NUMERIC_TYPE nodata_elevation,
	int acceleration, 
	int verbose
)
{
	_0 = cuda::Topography::load(filename, geometry, pitch, offset,
			nodata_elevation, acceleration, verbose); 
	_1x = load_slope(filename, "1x", geometry, pitch, offset, verbose);
	_1y = load_slope(filename, "1y", geometry, pitch, offset, verbose);

	zero_perimeter_slopes(geometry, pitch);
	clamp_boundary_values(geometry, pitch);
}

NUMERIC_TYPE* lis::cuda::dg2::Topography::load_slope
(
	const char* filename,
	const char* suffix,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose
)
{
	char dem_slope[256]; 
	strcpy(dem_slope, filename);
	strcat(dem_slope, suffix); 

	NUMERIC_TYPE* array = lis::GhostRaster::allocate(geometry);

	FILE* file = fopen_or_die(dem_slope, "rb", "Loading DEM slopes\n", verbose);
	Geometry slope_geometry;
	NUMERIC_TYPE no_data_value;
	AsciiRaster::read_header(file, slope_geometry, no_data_value);
	AsciiRaster::match_cell_dimensions_or_die(geometry, slope_geometry,
			"dg2::Topography::load");
	AsciiRaster::read(file, array, geometry, pitch, offset);
	fclose(file);

	return array;
};

void lis::cuda::dg2::Topography::zero_perimeter_slopes
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
			_1x[j*pitch + i + 1] = C(0.0);
			_1y[j * pitch + i + 1] = C(0.0);
		}

		// east
		{
			const int i=geometry.xsz+1;
			_1x[j*pitch + i - 1] = C(0.0);
			_1y[j * pitch + i - 1] = C(0.0);
		}
	}

	for (int i=1; i<geometry.xsz+1; i++)
	{
		// north
		{
			const int j=0;
			_1y[(j+1)*pitch + i] = C(0.0);
			_1x[(j + 1) * pitch + i] = C(0.0);
		}

		// south
		{
			const int j=geometry.ysz+1;
			_1y[(j-1)*pitch + i] = C(0.0);
			_1x[(j - 1) * pitch + i] = C(0.0);
		}
	}
}

void lis::cuda::dg2::Topography::clamp_boundary_values
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
			_0[j*pitch + i] = _0[j*pitch + i+1];
				//- SQRT(C(3.0))*_1x[j*pitch + i+1];
		}

		// east
		{
			const int i=geometry.xsz+1;
			_0[j * pitch + i] = _0[j * pitch + i - 1];
				//+ SQRT(C(3.0))*_1x[j*pitch + i-1];
		}
	}

	for (int i=1; i<geometry.xsz+1; i++)
	{
		// north
		{
			const int j=0;
			_0[j * pitch + i] = _0[(j + 1) * pitch + i];
				//+ SQRT(C(3.0))*_1y[(j+1)*pitch + i];
		}

		// south
		{
			const int j=geometry.ysz+1;
			_0[j * pitch + i] = _0[(j - 1) * pitch + i];
				//- SQRT(C(3.0))*_1y[(j-1)*pitch + i];
		}
	}
}

lis::cuda::dg2::Topography::~Topography()
{
	delete[] _0;
	delete[] _1x;
	delete[] _1y;
}

void lis::cuda::dg2::DeviceTopography::initialise
(
	Topography& DEM,
	Geometry& geometry
)
{
	_0 = cuda::GhostRaster::allocate_device(geometry);
	_1x = cuda::GhostRaster::allocate_device(geometry);
	_1y = cuda::GhostRaster::allocate_device(geometry);
	cuda::GhostRaster::copy(_0, DEM._0, geometry);
	cuda::GhostRaster::copy(_1x, DEM._1x, geometry);
	cuda::GhostRaster::copy(_1y, DEM._1y, geometry);

	Zstar_x = cuda::GhostRaster::allocate_device(geometry);
	Zstar_y = cuda::GhostRaster::allocate_device(geometry);
	initialise_Zstar_x(); 
	initialise_Zstar_y(); 
}

void lis::cuda::dg2::DeviceTopography::initialise_Zstar
(

)
{

	initialise_Zstar_x();
	initialise_Zstar_y();
}

__device__
NUMERIC_TYPE lis::cuda::dg2::DeviceTopography::neg_x
(
	int idx
)
{
	return _0[idx] + SQRT(C(3.0))*_1x[idx];
}

__device__
NUMERIC_TYPE lis::cuda::dg2::DeviceTopography::pos_x
(
	int idx
)
{
	return _0[idx] - SQRT(C(3.0))*_1x[idx];
}

__device__
NUMERIC_TYPE lis::cuda::dg2::DeviceTopography::neg_y
(
	int idx
)
{
	return _0[idx] + SQRT(C(3.0))*_1y[idx];
}

__device__
NUMERIC_TYPE lis::cuda::dg2::DeviceTopography::pos_y
(
	int idx
)
{
	return _0[idx] - SQRT(C(3.0))*_1y[idx];
}

__device__
NUMERIC_TYPE lis::cuda::dg2::DeviceTopography::Zdagger
(
	NUMERIC_TYPE ETA,
	NUMERIC_TYPE Zstar
)
{
	return Zstar - FMAX(C(0.0), -(ETA - Zstar));
}

void lis::cuda::dg2::DeviceTopography::initialise_Zstar_x()
{
	lis::cuda::dg2::initialise_Zstar_x<<<1, cuda::block_size>>>(Zstar_x, *this);
}

void lis::cuda::dg2::DeviceTopography::initialise_Zstar_y()
{
	lis::cuda::dg2::initialise_Zstar_y<<<1, cuda::block_size>>>(Zstar_y, *this);
}

void lis::cuda::dg2::DeviceTopography::free_device()
{
	cuda::free_device(_0);
	cuda::free_device(_1x);
	cuda::free_device(_1y);
}
