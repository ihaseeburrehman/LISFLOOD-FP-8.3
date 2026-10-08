#include "cuda_utils.cuh"

cudaError_t lis::cuda::acc_nugrid::sync()
{
	return cudaDeviceSynchronize();
}

cudaError_t lis::cuda::acc_nugrid::peek()
{
	return cudaPeekAtLastError();
}

cudaError_t lis::cuda::acc_nugrid::reset()
{
	return cudaDeviceReset();
}

cudaError_t lis::cuda::acc_nugrid::copy_cuda
(
	void* dst,
	void* src,
	size_t bytes
)
{
	cudaError_t error = cudaMemcpy
	(
		dst,
		src,
		bytes,
		cudaMemcpyDefault
	);

	return error;
}

__host__ __device__
void* lis::cuda::acc_nugrid::malloc_device
(
	size_t bytes
)
{
	void* ptr;
	
	cudaMalloc
	(
		&ptr,
		bytes
	);

	return ptr;
}

__host__ __device__
void* lis::cuda::acc_nugrid::malloc_unified
(
	size_t bytes
)
{
	void* ptr;

	cudaMallocManaged
	(
		&ptr,
		bytes
	);

	return ptr;
}


__host__ __device__
cudaError_t lis::cuda::acc_nugrid::free_device
(
	void* ptr
)
{
	return (nullptr != ptr) ? cudaFree(ptr) : cudaSuccess;
}

__host__ __device__
cudaError_t lis::cuda::acc_nugrid::free_unified
(
	void* ptr
)
{
	return (nullptr != ptr) ? cudaFree(ptr) : cudaSuccess;
}