#pragma once
#include "cuda_dem.cuh"

namespace lis
{
namespace cuda
{
namespace dg2
{

class Topography
{
public:
	Topography
	(
		const char* filename,
		Geometry& geometry,
		int& pitch,
		int& offset,
		NUMERIC_TYPE nodata_elevation,
		int acceleration, 
		int verbose
	);

	~Topography();

	NUMERIC_TYPE* __restrict__ _0;
	NUMERIC_TYPE* __restrict__ _1x;
	NUMERIC_TYPE* __restrict__ _1y;

private:
	NUMERIC_TYPE* load_slope
	(
		const char* filename,
		const char* suffix,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose
	);

	void zero_perimeter_slopes
	(
		Geometry& geometry,
		int pitch
	);

	void clamp_boundary_values
	(
		Geometry& geometry,
		int pitch
	);
};

class DeviceTopography
{
public:
	void initialise
	(
		Topography& DEM,
		Geometry& geometry
	);
	
	void initialise_Zstar
	(

	);

	__device__
	NUMERIC_TYPE neg_x
	(
		int idx
	);

	__device__
	NUMERIC_TYPE pos_x
	(
		int idx
	);

	__device__
	NUMERIC_TYPE neg_y
	(
		int idx
	);

	__device__
	NUMERIC_TYPE pos_y
	(
		int idx
	);

	__device__
	NUMERIC_TYPE Zdagger
	(
		NUMERIC_TYPE ETA,
		NUMERIC_TYPE Zstar
	);

	void free_device();

	NUMERIC_TYPE* __restrict__ _0;
	NUMERIC_TYPE* __restrict__ _1x;
	NUMERIC_TYPE* __restrict__ _1y;
	NUMERIC_TYPE* __restrict__ Zstar_x;
	NUMERIC_TYPE* __restrict__ Zstar_y;

private:
	void initialise_Zstar_x();
	void initialise_Zstar_y();
};

}
}
}
