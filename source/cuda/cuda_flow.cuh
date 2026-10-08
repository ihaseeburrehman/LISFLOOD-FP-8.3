#pragma once
#include "../lisflood.h"

namespace lis
{
namespace cuda
{

struct FlowVector
{
	NUMERIC_TYPE H;
	NUMERIC_TYPE HU;
	NUMERIC_TYPE HV;

	__device__ FlowVector star
	(
		NUMERIC_TYPE Z,
		NUMERIC_TYPE Zstar
	) const;

	__host__ __device__
	FlowVector operator-() const
	{
		return { -H, -HU, -HV };
	}

	__device__
	FlowVector physical_flux_x() const;

	__device__
	FlowVector physical_flux_y() const;

private:
	__device__ NUMERIC_TYPE calculate_Hstar
	(
		NUMERIC_TYPE Z,
		NUMERIC_TYPE Zstar
	) const;

	__device__ NUMERIC_TYPE speed
	(
		NUMERIC_TYPE discharge
	) const;
};

__host__ __device__
inline FlowVector operator+
(
	const FlowVector& lhs,
	const FlowVector& rhs
)
{
	return { lhs.H+rhs.H, lhs.HU+rhs.HU, lhs.HV+rhs.HV };
}

__host__ __device__
inline FlowVector operator-
(
	const FlowVector& lhs,
	const FlowVector& rhs
)
{
	return { lhs.H-rhs.H, lhs.HU-rhs.HU, lhs.HV-rhs.HV };
}

__host__ __device__
inline FlowVector operator*
(
	NUMERIC_TYPE lhs,
	const FlowVector& rhs
)
{
	return { lhs*rhs.H, lhs*rhs.HU, lhs*rhs.HV };
}

__host__ __device__
inline FlowVector operator/
(
	const FlowVector& lhs,
	NUMERIC_TYPE rhs
)
{
	return { lhs.H/rhs, lhs.HU/rhs, lhs.HV/rhs };
}

}
}
