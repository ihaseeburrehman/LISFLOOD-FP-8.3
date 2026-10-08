#include "ghostraster.h"
#include "../geometry.h"

int lis::GhostRaster::elements(Geometry& geometry)
{
	return (geometry.xsz+2)*(geometry.ysz+2);
}

int lis::GhostRaster::elements_H(Geometry& geometry)
{
	return (geometry.xsz)*(geometry.ysz);
}

int lis::GhostRaster::elements_Q(Geometry& geometry)
{
	return (geometry.xsz + 1)*(geometry.ysz + 1);
}

int lis::GhostRaster::pitch(Geometry& geometry)
{
	return geometry.xsz + 2;
}

int lis::GhostRaster::offset(Geometry& geometry)
{
	return geometry.xsz + 3;
}

int lis::GhostRaster::pitch_ACC(Geometry& geometry)
{
	return geometry.xsz;
}

int lis::GhostRaster::offset_ACC(Geometry& geometry)
{
	return 0;
}

NUMERIC_TYPE* lis::GhostRaster::allocate(Geometry& geometry)
{
	return new NUMERIC_TYPE[elements(geometry)]();
}

NUMERIC_TYPE* lis::GhostRaster::allocate_H(Geometry& geometry)
{
	return new NUMERIC_TYPE[elements_H(geometry)]();
}