#include "ghostraster.h"
#include "cuda_geometry.cuh"
#include "cuda_util.cuh"

NUMERIC_TYPE* lis::cuda::GhostRaster::allocate_pinned
(
	Geometry& geometry
)
{
	return static_cast<NUMERIC_TYPE*>(malloc_pinned(
				lis::GhostRaster::elements(geometry)*sizeof(NUMERIC_TYPE)));
}

NUMERIC_TYPE* lis::cuda::GhostRaster::allocate_device
(
	Geometry& geometry
)
{
	return static_cast<NUMERIC_TYPE*>(malloc_device(
				lis::GhostRaster::elements(geometry)*sizeof(NUMERIC_TYPE)));
}

void lis::cuda::GhostRaster::copy
(
	NUMERIC_TYPE* dst,
	NUMERIC_TYPE* src,
	Geometry& geometry
)
{
	checkCudaErrors(cudaMemcpy(dst, src,
				lis::GhostRaster::elements(geometry)*sizeof(NUMERIC_TYPE),
				cudaMemcpyDefault));
}



NUMERIC_TYPE* lis::cuda::GhostRaster::allocate_pinned_H
(
	Geometry& geometry
)
{
	return static_cast<NUMERIC_TYPE*>(malloc_pinned(
				lis::GhostRaster::elements_H(geometry)*sizeof(NUMERIC_TYPE)));
}

NUMERIC_TYPE* lis::cuda::GhostRaster::allocate_device_H
(
	Geometry& geometry
)
{
	return static_cast<NUMERIC_TYPE*>(malloc_device(
				lis::GhostRaster::elements_H(geometry)*sizeof(NUMERIC_TYPE)));
}

void lis::cuda::GhostRaster::copy_H
(
	NUMERIC_TYPE* dst,
	NUMERIC_TYPE* src,
	Geometry& geometry
)
{
	checkCudaErrors(cudaMemcpy(dst, src,
				lis::GhostRaster::elements_H(geometry)*sizeof(NUMERIC_TYPE),
				cudaMemcpyDefault));
}


NUMERIC_TYPE* lis::cuda::GhostRaster::allocate_pinned_Q
(
	Geometry& geometry
)
{
	return static_cast<NUMERIC_TYPE*>(malloc_pinned(
				lis::GhostRaster::elements_Q(geometry)*sizeof(NUMERIC_TYPE)));
}

NUMERIC_TYPE* lis::cuda::GhostRaster::allocate_device_Q
(
	Geometry& geometry
)
{
	return static_cast<NUMERIC_TYPE*>(malloc_device(
				lis::GhostRaster::elements_Q(geometry)*sizeof(NUMERIC_TYPE)));
}

void lis::cuda::GhostRaster::copy_Q
(
	NUMERIC_TYPE* dst,
	NUMERIC_TYPE* src,
	Geometry& geometry
)
{
	checkCudaErrors(cudaMemcpy(dst, src,
				lis::GhostRaster::elements_Q(geometry)*sizeof(NUMERIC_TYPE),
				cudaMemcpyDefault));
}