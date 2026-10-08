#include "cuda_snapshot.cuh"
#include "io.h"

lis::cuda::VelocityWriter::VelocityWriter
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* discharge,
	NUMERIC_TYPE DepthThresh
)
:
H(H),
discharge(discharge),
DepthThresh(DepthThresh)
{}

NUMERIC_TYPE lis::cuda::VelocityWriter::operator[]
(
	int idx
)
{
	if (H[idx] > DepthThresh)
	{
		return discharge[idx] / H[idx];
	}
	else
	{
		return C(0.0);
	}
}

lis::cuda::ElevationWriter::ElevationWriter
(
	const NUMERIC_TYPE* H,
	const NUMERIC_TYPE* DEM
)
:
H(H),
DEM(DEM)
{}

NUMERIC_TYPE lis::cuda::ElevationWriter::operator[]
(
	int idx
) const
{
	return DEM[idx] + H[idx];
}
