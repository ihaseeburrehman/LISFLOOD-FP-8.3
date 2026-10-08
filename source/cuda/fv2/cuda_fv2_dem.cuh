#pragma once
#include "cuda_dem.cuh"

namespace lis
{
namespace cuda
{
namespace fv2
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
		int acceleration, //// added
		int verbose
	);

	~Topography();

	NUMERIC_TYPE* __restrict__ _0;

private:

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

	__device__
	NUMERIC_TYPE Zdagger
	(
		NUMERIC_TYPE ETA,
		NUMERIC_TYPE Zstar
	);

	void free_device();

	NUMERIC_TYPE* __restrict__ _0;
	NUMERIC_TYPE* __restrict__ Zstar_x;
	NUMERIC_TYPE* __restrict__ Zstar_y;

private:
	void initialise_Zstar_x();
	void initialise_Zstar_y();
};

}
}
}
