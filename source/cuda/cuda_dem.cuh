#pragma once
#include "../lisflood.h"
#include "../geometry.h"

namespace lis
{
namespace cuda
{

struct Topography
{
	static NUMERIC_TYPE* load
	(
		const char* filename,
		Geometry& geometry,
		int& pitch,
		int& offset,
		NUMERIC_TYPE nodata_elevation,
		int acceleration,
		int verbose
	);

	static void initialise_Zstar_x
	(
		NUMERIC_TYPE* __restrict__ Zstar_x,
		const NUMERIC_TYPE* __restrict__ DEM
	);

	static void initialise_Zstar_y
	(
		NUMERIC_TYPE* __restrict__ Zstar_y,
		const NUMERIC_TYPE* __restrict__ DEM
	);

	static void clamp_boundary_values
	(
		NUMERIC_TYPE* DEM,
		Geometry& geometry,
		int pitch
	);
};

}
}
