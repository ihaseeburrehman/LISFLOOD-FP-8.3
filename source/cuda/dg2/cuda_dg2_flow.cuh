#pragma once
#include "../../lisflood.h"
#include "cuda_flow.cuh"
#include "cuda_geometry.cuh"

namespace lis
{
namespace cuda
{
namespace dg2
{

struct FlowCoeffs;
struct FlowCoeffsRef;

struct FlowCoeffs
{
	NUMERIC_TYPE H;
	NUMERIC_TYPE H1x;
	NUMERIC_TYPE H1y;
	NUMERIC_TYPE HU;
	NUMERIC_TYPE HU1x;
	NUMERIC_TYPE HU1y;
	NUMERIC_TYPE HV;
	NUMERIC_TYPE HV1x;
	NUMERIC_TYPE HV1y;

	FlowCoeffs() = default;

	__host__ __device__
	FlowCoeffs
	(
		NUMERIC_TYPE H,
		NUMERIC_TYPE H1x,
		NUMERIC_TYPE H1y,
		NUMERIC_TYPE HU,
		NUMERIC_TYPE HU1x,
		NUMERIC_TYPE HU1y,
		NUMERIC_TYPE HV,
		NUMERIC_TYPE HV1x,
		NUMERIC_TYPE HV1y
	);

	__host__ __device__
	FlowCoeffs
	(
		const FlowCoeffsRef& src
	);

	__device__
	void set_0
	(
		const FlowVector& v
	);

	__device__
	void set_1x
	(
		const FlowVector& v
	);

	__device__
	void set_1y
	(
		const FlowVector& v
	);

	__device__
	FlowVector get_0();

	__device__
	FlowVector neg_x();

	__device__
	FlowVector pos_x();

	__device__
	FlowVector neg_y();

	__device__
	FlowVector pos_y();

	__device__
	FlowVector gauss_lower_x();

	__device__
	FlowVector gauss_upper_x();

	__device__
	FlowVector gauss_lower_y();

	__device__
	FlowVector gauss_upper_y();

	__device__
	FlowCoeffs& operator+=
	(
		const FlowCoeffs& rhs
	);
};

__host__ __device__
inline FlowCoeffs operator+
(
	const FlowCoeffs& lhs,
	const FlowCoeffs& rhs
)
{
	return { lhs.H+rhs.H, lhs.H1x+rhs.H1x, lhs.H1y+rhs.H1y,
		lhs.HU+rhs.HU, lhs.HU1x+rhs.HU1x, lhs.HU1y+rhs.HU1y,
		lhs.HV+rhs.HV, lhs.HV1x+rhs.HV1x, lhs.HV1y+rhs.HV1y };
}

__host__ __device__
inline FlowCoeffs operator*
(
	NUMERIC_TYPE lhs,
	const FlowCoeffs& rhs
)
{
	return { lhs*rhs.H, lhs*rhs.H1x, lhs*rhs.H1y,
		lhs*rhs.HU, lhs*rhs.HU1x, lhs*rhs.HU1y,
		lhs*rhs.HV, lhs*rhs.HV1x, lhs*rhs.HV1y };
}

struct FlowCoeffsRef
{
	NUMERIC_TYPE& H;
	NUMERIC_TYPE& H1x;
	NUMERIC_TYPE& H1y;
	NUMERIC_TYPE& HU;
	NUMERIC_TYPE& HU1x;
	NUMERIC_TYPE& HU1y;
	NUMERIC_TYPE& HV;
	NUMERIC_TYPE& HV1x;
	NUMERIC_TYPE& HV1y;

	__host__ __device__
	FlowCoeffsRef& operator=
	(
		const FlowCoeffs& rhs
	);

	__device__
	void set_0
	(
		const FlowVector& v
	);
};

struct Flow
{
	NUMERIC_TYPE* __restrict__ H;
	NUMERIC_TYPE* __restrict__ H1x;
	NUMERIC_TYPE* __restrict__ H1y;
	NUMERIC_TYPE* __restrict__ HU;
	NUMERIC_TYPE* __restrict__ HU1x;
	NUMERIC_TYPE* __restrict__ HU1y;
	NUMERIC_TYPE* __restrict__ HV;
	NUMERIC_TYPE* __restrict__ HV1x;
	NUMERIC_TYPE* __restrict__ HV1y;

	NUMERIC_TYPE* __restrict__ maxH;
	NUMERIC_TYPE* __restrict__ totalHtm;
	NUMERIC_TYPE* __restrict__ maxHtm;
	NUMERIC_TYPE* __restrict__ initHtm;

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

	__host__ __device__ FlowCoeffsRef operator[]
	(
		int idx
	)
	{
		return { H[idx], H1x[idx], H1y[idx], HU[idx], HU1x[idx], HU1y[idx],
			HV[idx], HV1x[idx], HV1y[idx] };
	}

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
