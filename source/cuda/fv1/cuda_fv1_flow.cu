#include "cuda_fv1_flow.cuh"
#include "cuda_geometry.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"

void lis::cuda::fv1::Flow::allocate_pinned
(
	Flow& flow,
	Geometry& geometry
)
{
	flow.H = cuda::GhostRaster::allocate_pinned(geometry);
	flow.HU = cuda::GhostRaster::allocate_pinned(geometry);
	flow.HV = cuda::GhostRaster::allocate_pinned(geometry);
}

void lis::cuda::fv1::Flow::allocate_device
(
	Flow& flow,
	Geometry& geometry
)
{
	flow.H = cuda::GhostRaster::allocate_device(geometry);
	flow.HU = cuda::GhostRaster::allocate_device(geometry);
	flow.HV = cuda::GhostRaster::allocate_device(geometry);
}

void lis::cuda::fv1::Flow::copy
(
	Flow& dst,
	Flow& src,
	Geometry& geometry
)
{
	cuda::GhostRaster::copy(dst.H, src.H, geometry);
	cuda::GhostRaster::copy(dst.HU, src.HU, geometry);
	cuda::GhostRaster::copy(dst.HV, src.HV, geometry);
}

void lis::cuda::fv1::Flow::free_pinned
(
	Flow& flow
)
{
	cuda::free_pinned(flow.H);
	cuda::free_pinned(flow.HU);
	cuda::free_pinned(flow.HV);
}

void lis::cuda::fv1::Flow::free_device
(
	Flow& flow
)
{
	cuda::free_device(flow.H);
	cuda::free_device(flow.HU);
	cuda::free_device(flow.HV);
}
