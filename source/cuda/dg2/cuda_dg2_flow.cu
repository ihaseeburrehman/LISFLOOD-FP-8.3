#include "cuda_dg2_flow.cuh"
#include "cuda_util.cuh"

__host__ __device__
lis::cuda::dg2::FlowCoeffs::FlowCoeffs
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
lis::cuda::dg2::FlowCoeffs::FlowCoeffs
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
void lis::cuda::dg2::FlowCoeffs::set_0
(
	const FlowVector& v
)
{
	H = v.H;
	HU = v.HU;
	HV = v.HV;
}

__device__
void lis::cuda::dg2::FlowCoeffs::set_1x
(
	const FlowVector& v
)
{
	H1x = v.H;
	HU1x = v.HU;
	HV1x = v.HV;
}

__device__
void lis::cuda::dg2::FlowCoeffs::set_1y
(
	const FlowVector& v
)
{
	H1y = v.H;
	HU1y = v.HU;
	HV1y = v.HV;
}

__device__
lis::cuda::FlowVector lis::cuda::dg2::FlowCoeffs::get_0()
{
	return { H, HU, HV };
}

__device__
lis::cuda::FlowVector lis::cuda::dg2::FlowCoeffs::neg_x()
{
	return { H + SQRT(C(3.0))*H1x, HU + SQRT(C(3.0))*HU1x,
		HV + SQRT(C(3.0))*HV1x };
}

__device__
lis::cuda::FlowVector lis::cuda::dg2::FlowCoeffs::pos_x()
{
	return { H - SQRT(C(3.0))*H1x, HU - SQRT(C(3.0))*HU1x,
		HV - SQRT(C(3.0))*HV1x };
}

__device__
lis::cuda::FlowVector lis::cuda::dg2::FlowCoeffs::neg_y()
{
	return { H + SQRT(C(3.0))*H1y, HU + SQRT(C(3.0))*HU1y,
		HV + SQRT(C(3.0))*HV1y };
}

__device__
lis::cuda::FlowVector lis::cuda::dg2::FlowCoeffs::pos_y()
{
	return { H - SQRT(C(3.0))*H1y, HU - SQRT(C(3.0))*HU1y,
		HV - SQRT(C(3.0))*HV1y };
}

__device__
lis::cuda::FlowVector lis::cuda::dg2::FlowCoeffs::gauss_lower_x()
{
	return { H - H1x, HU - HU1x, HV - HV1x };
}

__device__
lis::cuda::FlowVector lis::cuda::dg2::FlowCoeffs::gauss_upper_x()
{
	return { H + H1x, HU + HU1x, HV + HV1x };
}

__device__
lis::cuda::FlowVector lis::cuda::dg2::FlowCoeffs::gauss_lower_y()
{
	return { H - H1y, HU - HU1y, HV - HV1y };
}

__device__
lis::cuda::FlowVector lis::cuda::dg2::FlowCoeffs::gauss_upper_y()
{
	return { H + H1y, HU + HU1y, HV + HV1y };
}

__device__
lis::cuda::dg2::FlowCoeffs& lis::cuda::dg2::FlowCoeffs::operator+=
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
lis::cuda::dg2::FlowCoeffsRef& lis::cuda::dg2::FlowCoeffsRef::operator=
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
void lis::cuda::dg2::FlowCoeffsRef::set_0
(
	const FlowVector& v
)
{
	H = v.H;
	HU = v.HU;
	HV = v.HV;
}

void lis::cuda::dg2::Flow::allocate_pinned
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

void lis::cuda::dg2::Flow::allocate_device
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

void lis::cuda::dg2::Flow::copy
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

void lis::cuda::dg2::Flow::free_pinned
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

void lis::cuda::dg2::Flow::free_device
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

