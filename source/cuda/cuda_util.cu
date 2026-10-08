#include "cuda_util.cuh"
#include <helper_cuda.h>

void lis::cuda::sync()
{
	checkCudaErrors(cudaDeviceSynchronize());
}


void lis::cuda::peek()
{
	checkCudaErrors(cudaPeekAtLastError());
}


//cudaError_t lis::cuda::sync()  // for debug
//{
//	return cudaDeviceSynchronize();
//}
//
//cudaError_t lis::cuda::peek()
//{
//	return cudaPeekAtLastError();
//}


void lis::cuda::copy
(
	void* dst,
	void* src,
	size_t count
)
{
	checkCudaErrors(cudaMemcpy(dst, src, count, cudaMemcpyDefault));
}

void* lis::cuda::malloc_unified
(
	size_t size
)
{
	void* ptr;
	checkCudaErrors(cudaMallocManaged(&ptr, size));
	return ptr;
}

void* lis::cuda::malloc_pinned
(
	size_t size
)
{
	void* ptr;
	checkCudaErrors(cudaMallocHost(&ptr, size));
	return ptr;
}

void* lis::cuda::malloc_device
(
	size_t size
)
{
	void* ptr;
	checkCudaErrors(cudaMalloc(&ptr, size));
	return ptr;
}

void lis::cuda::free_unified
(
	void* ptr
)
{
	checkCudaErrors(cudaFree(ptr));
}

void lis::cuda::free_pinned
(
	void* ptr
)
{
	checkCudaErrors(cudaFreeHost(ptr));
}

void lis::cuda::free_device
(
	void* ptr
)
{
	checkCudaErrors(cudaFree(ptr));
}

int lis::cuda::get_device()
{
	int device;
	checkCudaErrors(cudaGetDevice(&device));
	return device;
}

void lis::cuda::get_device_properties
(
	cudaDeviceProp& properties,
	int device
)
{
	checkCudaErrors(cudaGetDeviceProperties(&properties, device));
}
