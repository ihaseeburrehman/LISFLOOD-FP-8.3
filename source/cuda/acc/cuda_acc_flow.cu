#include "cuda_acc_flow.cuh"
#include "cuda_geometry.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"

void lis::cuda::acc::Flow::allocate_pinned
(
	Flow& flow,
	Geometry& geometry
)
{
	flow.H = cuda::GhostRaster::allocate_pinned_H(geometry);
//	flow.ChanMask = cuda::GhostRaster::allocate_pinned_H(geometry);

	flow.maxH = cuda::GhostRaster::allocate_pinned_H(geometry);
	flow.totalHtm = cuda::GhostRaster::allocate_pinned_H(geometry);
	flow.maxHtm = cuda::GhostRaster::allocate_pinned_H(geometry);
	flow.initHtm = cuda::GhostRaster::allocate_pinned_H(geometry);

	flow.Qx = cuda::GhostRaster::allocate_pinned_Q(geometry);
	flow.Qy = cuda::GhostRaster::allocate_pinned_Q(geometry);
	flow.Qxold = cuda::GhostRaster::allocate_pinned_Q(geometry);
	flow.Qyold = cuda::GhostRaster::allocate_pinned_Q(geometry);
	flow.Vx = cuda::GhostRaster::allocate_pinned_Q(geometry);
	flow.Vy = cuda::GhostRaster::allocate_pinned_Q(geometry);
}

void lis::cuda::acc::Flow::allocate_device
(
	Flow& flow,
	Geometry& geometry
)
{
	flow.H = cuda::GhostRaster::allocate_device_H(geometry);
//	flow.ChanMask = cuda::GhostRaster::allocate_device_H(geometry);

	flow.maxH = cuda::GhostRaster::allocate_device_H(geometry);
	flow.totalHtm = cuda::GhostRaster::allocate_device_H(geometry);
	flow.maxHtm = cuda::GhostRaster::allocate_device_H(geometry);
	flow.initHtm = cuda::GhostRaster::allocate_device_H(geometry);

	
	flow.Qx = cuda::GhostRaster::allocate_device_Q(geometry);
	flow.Qy = cuda::GhostRaster::allocate_device_Q(geometry);
	flow.Qxold = cuda::GhostRaster::allocate_device_Q(geometry);
	flow.Qyold = cuda::GhostRaster::allocate_device_Q(geometry);
	flow.Vx = cuda::GhostRaster::allocate_device_Q(geometry);
	flow.Vy = cuda::GhostRaster::allocate_device_Q(geometry);
}

void lis::cuda::acc::Flow::copy
(
	Flow& dst,
	Flow& src,
	Geometry& geometry
)
{
	cuda::GhostRaster::copy_H(dst.H, src.H, geometry);
//	cuda::GhostRaster::copy_H(dst.ChanMask, src.ChanMask, geometry);

	cuda::GhostRaster::copy_H(dst.maxH, src.maxH, geometry);
	cuda::GhostRaster::copy_H(dst.totalHtm, src.totalHtm, geometry);
	cuda::GhostRaster::copy_H(dst.maxHtm, src.maxHtm, geometry);
	cuda::GhostRaster::copy_H(dst.initHtm, src.initHtm, geometry);

	
	cuda::GhostRaster::copy_Q(dst.Qx, src.Qx, geometry);
	cuda::GhostRaster::copy_Q(dst.Qy, src.Qy, geometry);
	cuda::GhostRaster::copy_Q(dst.Qxold, src.Qxold, geometry);
	cuda::GhostRaster::copy_Q(dst.Qyold, src.Qyold, geometry);
	cuda::GhostRaster::copy_Q(dst.Vx, src.Vx, geometry);
	cuda::GhostRaster::copy_Q(dst.Vy, src.Vy, geometry);
}

void lis::cuda::acc::Flow::free_pinned
(
	Flow& flow
)
{
	cuda::free_pinned(flow.H);
//	cuda::free_pinned(flow.ChanMask);

	cuda::free_pinned(flow.maxH);
	cuda::free_pinned(flow.totalHtm);
	cuda::free_pinned(flow.maxHtm);
	cuda::free_pinned(flow.initHtm);
	
	cuda::free_pinned(flow.Qx);
	cuda::free_pinned(flow.Qy);
	cuda::free_pinned(flow.Qxold);
	cuda::free_pinned(flow.Qyold);
	cuda::free_pinned(flow.Vx);
	cuda::free_pinned(flow.Vy);
}

void lis::cuda::acc::Flow::free_device
(
	Flow& flow
)
{
	cuda::free_device(flow.H); 
//	cuda::free_device(flow.ChanMask);

	cuda::free_device(flow.maxH);
	cuda::free_device(flow.totalHtm);
	cuda::free_device(flow.maxHtm);
	cuda::free_device(flow.initHtm);
	
	cuda::free_device(flow.Qx);
	cuda::free_device(flow.Qy);
	cuda::free_device(flow.Qxold);
	cuda::free_device(flow.Qyold);
	cuda::free_device(flow.Vx);
	cuda::free_device(flow.Vy);
}
