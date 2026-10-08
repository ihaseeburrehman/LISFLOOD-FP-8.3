#pragma once
#include "../geometry.h"
#include "cuda_flow.cuh"

namespace lis
{
namespace cuda
{
namespace fv1
{

struct Flow
{
	NUMERIC_TYPE* __restrict__ H;
	NUMERIC_TYPE* __restrict__ HU;
	NUMERIC_TYPE* __restrict__ HV;

	NUMERIC_TYPE* __restrict__ maxH;
	NUMERIC_TYPE* __restrict__ totalHtm;
	NUMERIC_TYPE* __restrict__ maxHtm;
	NUMERIC_TYPE* __restrict__ initHtm;

	__host__ __device__ FlowVector operator[]
	(
		int idx
	)
	{
		return { H[idx], HU[idx], HV[idx] };
	}

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
