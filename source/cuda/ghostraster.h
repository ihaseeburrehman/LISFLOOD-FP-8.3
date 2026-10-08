#pragma once
#include "../lisflood.h"
#include "../geometry.h"

namespace lis
{
	
struct GhostRaster
{
	static int elements(Geometry& geometry);
	static int elements_H(Geometry& geometry);
	static int elements_Q(Geometry& geometry);
	static int pitch(Geometry& geometry);
	static int offset(Geometry& geometry);
	static int pitch_ACC(Geometry& geometry);
	static int offset_ACC(Geometry& geometry);	
	static NUMERIC_TYPE* allocate(Geometry& geometry);
	static NUMERIC_TYPE* allocate_H(Geometry& geometry);
};

}

