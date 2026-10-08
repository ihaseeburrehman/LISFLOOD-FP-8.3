#include "cuda_dem.cuh"
#include "cuda_solver.cuh"
#include "io.h"

namespace lis
{
namespace cuda
{
namespace fv1
{

__global__ void initialise_Zstar_x
(
	NUMERIC_TYPE* __restrict__ Zstar_x,
	const NUMERIC_TYPE* __restrict__ DEM
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j+1; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z_neg = DEM[j*cuda::pitch + i];
			NUMERIC_TYPE Z_pos = DEM[j*cuda::pitch + i+1];

			Zstar_x[j*cuda::pitch + i] = FMAX(Z_neg, Z_pos);
		}
	}
}

__global__ void initialise_Zstar_y
(
	NUMERIC_TYPE* __restrict__ Zstar_y,
	const NUMERIC_TYPE* __restrict__ DEM
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i+1; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z_neg = DEM[(j+1)*cuda::pitch + i];
			NUMERIC_TYPE Z_pos = DEM[j*cuda::pitch + i];

			Zstar_y[j*cuda::pitch + i] = FMAX(Z_neg, Z_pos);
		}
	}
}

}
}
}

NUMERIC_TYPE* lis::cuda::Topography::load
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
	FILE* dem_file = fopen_or_die(filename, "rb", "Loading DEM\n", verbose);
	NUMERIC_TYPE original_no_data_value = C(0.0);
	AsciiRaster::read_header(dem_file, geometry, original_no_data_value);
	NUMERIC_TYPE* DEM; 

    if (acceleration == 0){ 
		DEM = lis::GhostRaster::allocate(geometry); 
	    pitch = lis::GhostRaster::pitch(geometry);
	    offset = lis::GhostRaster::offset(geometry);
	} else {                                            
		DEM = lis::GhostRaster::allocate_H(geometry); 
	    pitch = lis::GhostRaster::pitch_ACC(geometry);   
	    offset = lis::GhostRaster::offset_ACC(geometry);	
	}
	AsciiRaster::read(dem_file, DEM, geometry, pitch, offset); 
	AsciiRaster::replace_no_data(DEM, geometry, pitch, offset,
			original_no_data_value, nodata_elevation);

	fclose(dem_file);


	return DEM;
}

void lis::cuda::Topography::initialise_Zstar_x
(
	NUMERIC_TYPE* __restrict__ Zstar_x,
	const NUMERIC_TYPE* __restrict__ DEM
)
{
	lis::cuda::fv1::initialise_Zstar_x<<<1, cuda::block_size>>>(Zstar_x, DEM);
}

void lis::cuda::Topography::initialise_Zstar_y
(
	NUMERIC_TYPE* __restrict__ Zstar_y,
	const NUMERIC_TYPE* __restrict__ DEM
)
{
	lis::cuda::fv1::initialise_Zstar_y<<<1, cuda::block_size>>>(Zstar_y, DEM);
}

void lis::cuda::Topography::clamp_boundary_values
(
	NUMERIC_TYPE* DEM,
	Geometry& geometry,
	int pitch
)
{
	for (int j=1; j<geometry.ysz+1; j++)
	{
		// west
		{
			const int i=0;
			DEM[j*pitch + i] = DEM[j*pitch + i+1];
		}

		// east
		{
			const int i=geometry.xsz+1;
			//if (j == 161 || j == 162 || j == 163 || j == 164) {
			//	DEM[j * pitch + i] = -C(20.5);
			//}
			//else
			//{
				DEM[j * pitch + i] = DEM[j * pitch + i - 1];
//			}
		}
	}

	for (int i=1; i<geometry.xsz+1; i++)
	{
		// north
		{
			const int j=0;
			DEM[j*pitch + i] = DEM[(j+1)*pitch + i];
		}

		// south
		{
			const int j=geometry.ysz+1;
			DEM[j*pitch + i] = DEM[(j-1)*pitch + i];
		}
	}
}
