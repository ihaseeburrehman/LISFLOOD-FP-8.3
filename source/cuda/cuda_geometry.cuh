#pragma once
#include "../lisflood.h"
#include "../geometry.h"

namespace lis
{
namespace cuda
{

struct GhostRaster
{
	static NUMERIC_TYPE* allocate_pinned
	(
		Geometry& geometry
	);

	static NUMERIC_TYPE* allocate_device
	(
		Geometry& geometry
	);

	static void copy
	(
		NUMERIC_TYPE* dst,
		NUMERIC_TYPE* src,
		Geometry& geometry
	);
	
	static NUMERIC_TYPE* allocate_pinned_H 
	(
		Geometry& geometry
	);

	static NUMERIC_TYPE* allocate_device_H 
	(
		Geometry& geometry
	);

	static void copy_H 
	(
		NUMERIC_TYPE* dst,
		NUMERIC_TYPE* src,
		Geometry& geometry
	);	
	
	static NUMERIC_TYPE* allocate_pinned_Q 
	(
		Geometry& geometry
	);

	static NUMERIC_TYPE* allocate_device_Q 
	(
		Geometry& geometry
	);

	static void copy_Q 
	(
		NUMERIC_TYPE* dst,
		NUMERIC_TYPE* src,
		Geometry& geometry
	);
	
};

}
}
