#pragma once

//#include "cuda_runtime.h" // for debug

namespace lis
{
namespace cuda
{

void sync();
void peek();

//cudaError_t sync();  // for debug
// 
//cudaError_t peek(); // for debug

void copy
(
	void* dst,
	void* src,
	size_t count
);

template<typename T>
void copy_to_symbol
(
	const T& symbol,
	const void* src,
	size_t count
);

void* malloc_unified
(
	size_t size
);

void* malloc_pinned
(
	size_t size
);

void* malloc_device
(
	size_t size
);

void free_unified
(
	void* ptr
);

void free_pinned
(
	void* ptr
);

void free_device
(
	void* ptr
);

int get_device();

void get_device_properties
(
	cudaDeviceProp& properties,
	int device
);

}
}

#include "cuda_util.templates.cu"
