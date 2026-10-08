#include "cuda_fv2_flow.cuh"
#include "cuda_util.cuh"

__host__ __device__
lis::cuda::fv2::FlowCoeffs::FlowCoeffs
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
)
:
	H(H),
	H1x(H1x),
	H1y(H1y),
	HU(HU),
	HU1x(HU1x),
	HU1y(HU1y),
	HV(HV),
	HV1x(HV1x),
	HV1y(HV1y)
{}

__host__ __device__
lis::cuda::fv2::FlowCoeffs::FlowCoeffs
(
	const FlowCoeffsRef& src
)
{
	H = src.H;
	H1x = src.H1x;
	H1y = src.H1y;
	HU = src.HU;
	HU1x = src.HU1x;
	HU1y = src.HU1y;
	HV = src.HV;
	HV1x = src.HV1x;
	HV1y = src.HV1y;
}

__device__
void lis::cuda::fv2::FlowCoeffs::set_0
(
	const FlowVector& v
)
{
	H = v.H;
	HU = v.HU;
	HV = v.HV;
}


__device__
lis::cuda::FlowVector lis::cuda::fv2::FlowCoeffs::get_0()
{
	return { H, HU, HV };
}

__device__
lis::cuda::FlowVector lis::cuda::fv2::FlowCoeffs::neg_x()
{
	return { H + C(0.5) * H1x, HU + C(0.5) * HU1x,
		HV + C(0.5) * HV1x };
}

__device__
lis::cuda::FlowVector lis::cuda::fv2::FlowCoeffs::pos_x()
{
	return { H - C(0.5) * H1x, HU - C(0.5) * HU1x,
		HV - C(0.5) * HV1x };
}

__device__
lis::cuda::FlowVector lis::cuda::fv2::FlowCoeffs::neg_y()
{
	return { H + C(0.5) * H1y, HU + C(0.5) * HU1y,
		HV + C(0.5) * HV1y };
}

__device__
lis::cuda::FlowVector lis::cuda::fv2::FlowCoeffs::pos_y()
{
	return { H - C(0.5) * H1y, HU - C(0.5) * HU1y,
		HV - C(0.5) * HV1y };
}


__device__
lis::cuda::fv2::FlowCoeffs& lis::cuda::fv2::FlowCoeffs::operator+=
(
	const FlowCoeffs& rhs
)
{
	H += rhs.H;
	H1x += rhs.H1x;
	H1y += rhs.H1y;
	HU += rhs.HU;
	HU1x += rhs.HU1x;
	HU1y += rhs.HU1y;
	HV += rhs.HV;
	HV1x += rhs.HV1x;
	HV1y += rhs.HV1y;

	return *this;
}

__host__ __device__
lis::cuda::fv2::FlowCoeffsRef& lis::cuda::fv2::FlowCoeffsRef::operator=
(
	const FlowCoeffs& rhs
)
{
	H = rhs.H;
	H1x = rhs.H1x;
	H1y = rhs.H1y;
	HU = rhs.HU;
	HU1x = rhs.HU1x;
	HU1y = rhs.HU1y;
	HV = rhs.HV;
	HV1x = rhs.HV1x;
	HV1y = rhs.HV1y;
	
	return *this;
}

__device__
void lis::cuda::fv2::FlowCoeffsRef::set_0
(
	const FlowVector& v
)
{
	H = v.H;
	HU = v.HU;
	HV = v.HV;
}

void lis::cuda::fv2::Flow::allocate_pinned
(
	Flow& flow,
	Geometry& geometry
)
{
	flow.H = cuda::GhostRaster::allocate_pinned(geometry);
	flow.H1x = cuda::GhostRaster::allocate_pinned(geometry);
	flow.H1y = cuda::GhostRaster::allocate_pinned(geometry);
	flow.HU = cuda::GhostRaster::allocate_pinned(geometry);
	flow.HU1x = cuda::GhostRaster::allocate_pinned(geometry);
	flow.HU1y = cuda::GhostRaster::allocate_pinned(geometry);
	flow.HV = cuda::GhostRaster::allocate_pinned(geometry);
	flow.HV1x = cuda::GhostRaster::allocate_pinned(geometry);
	flow.HV1y = cuda::GhostRaster::allocate_pinned(geometry);
}

void lis::cuda::fv2::Flow::allocate_device
(
	Flow& flow,
	Geometry& geometry
)
{
	flow.H = cuda::GhostRaster::allocate_device(geometry);
	flow.H1x = cuda::GhostRaster::allocate_device(geometry);
	flow.H1y = cuda::GhostRaster::allocate_device(geometry);
	flow.HU = cuda::GhostRaster::allocate_device(geometry);
	flow.HU1x = cuda::GhostRaster::allocate_device(geometry);
	flow.HU1y = cuda::GhostRaster::allocate_device(geometry);
	flow.HV = cuda::GhostRaster::allocate_device(geometry);
	flow.HV1x = cuda::GhostRaster::allocate_device(geometry);
	flow.HV1y = cuda::GhostRaster::allocate_device(geometry);
}

void lis::cuda::fv2::Flow::copy
(
	Flow& dst,
	Flow& src,
	Geometry& geometry
)
{
	cuda::GhostRaster::copy(dst.H, src.H, geometry);
	cuda::GhostRaster::copy(dst.H1x, src.H1x, geometry);
	cuda::GhostRaster::copy(dst.H1y, src.H1y, geometry);
	cuda::GhostRaster::copy(dst.HU, src.HU, geometry);
	cuda::GhostRaster::copy(dst.HU1x, src.HU1x, geometry);
	cuda::GhostRaster::copy(dst.HU1y, src.HU1y, geometry);
	cuda::GhostRaster::copy(dst.HV, src.HV, geometry);
	cuda::GhostRaster::copy(dst.HV1x, src.HV1x, geometry);
	cuda::GhostRaster::copy(dst.HV1y, src.HV1y, geometry);
}

void lis::cuda::fv2::Flow::free_pinned
(
	Flow& flow
)
{
	cuda::free_pinned(flow.H);
	cuda::free_pinned(flow.H1x);
	cuda::free_pinned(flow.H1y);
	cuda::free_pinned(flow.HU);
	cuda::free_pinned(flow.HU1x);
	cuda::free_pinned(flow.HU1y);
	cuda::free_pinned(flow.HV);
	cuda::free_pinned(flow.HV1x);
	cuda::free_pinned(flow.HV1y);
}

void lis::cuda::fv2::Flow::free_device
(
	Flow& flow
)
{
	cuda::free_device(flow.H);
	cuda::free_device(flow.H1x);
	cuda::free_device(flow.H1y);
	cuda::free_device(flow.HU);
	cuda::free_device(flow.HU1x);
	cuda::free_device(flow.HU1y);
	cuda::free_device(flow.HV);
	cuda::free_device(flow.HV1x);
	cuda::free_device(flow.HV1y);
}

