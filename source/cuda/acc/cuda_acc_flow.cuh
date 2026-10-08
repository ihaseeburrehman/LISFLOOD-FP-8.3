#pragma once
#include "../geometry.h"
#include "cuda_flow.cuh"

namespace lis
{
namespace cuda
{
namespace acc
{

struct Flow
{
	NUMERIC_TYPE* __restrict__ H;
//	NUMERIC_TYPE* __restrict__ ChanMask;

	NUMERIC_TYPE* __restrict__ maxH;
	NUMERIC_TYPE* __restrict__ totalHtm;
	NUMERIC_TYPE* __restrict__ maxHtm;
	NUMERIC_TYPE* __restrict__ initHtm;

	
	NUMERIC_TYPE* __restrict__ Qx;
	NUMERIC_TYPE* __restrict__ Qy;
	NUMERIC_TYPE* __restrict__ Qxold;
	NUMERIC_TYPE* __restrict__ Qyold;
	
	NUMERIC_TYPE* __restrict__ Vx;
	NUMERIC_TYPE* __restrict__ Vy;

	static void allocate_pinned
	(
		Flow& flow,
		Geometry& geometry
	);

	static void allocate_device
	(
		Flow& flow,
		Geometry& geometry
	);

	static void copy
	(
		Flow& dst,
		Flow& src,
		Geometry& geometry
	);

	static void free_pinned
	(
		Flow& flow
	);

	static void free_device
	(
		Flow& flow
	);
};
	
}
}
}
