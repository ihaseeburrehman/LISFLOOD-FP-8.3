#include "cuda_acc_max.cuh"
#include "cuda_acc_solver.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"
#include <helper_cuda.h>
#include <cub/cub.cuh>

lis::cuda::acc::MaxH::MaxH
(
	Geometry& geometry
)
:
elements(lis::GhostRaster::elements_H(geometry))
{
	temp = nullptr;
	NUMERIC_TYPE* dummy_in = nullptr;
	NUMERIC_TYPE* dummy_out = nullptr;

	checkCudaErrors(cub::DeviceReduce::Max(temp, bytes,
				dummy_in, dummy_out, elements));

	temp = malloc_device(bytes);
}

void lis::cuda::acc::MaxH::update
(
	const Flow& flow
)
{
	checkCudaErrors(cub::DeviceReduce::Max(temp, bytes, flow.H,
				&(cuda::acc::maxH), elements));
}

lis::cuda::acc::MaxH::~MaxH()
{
	free_device(temp);
}
