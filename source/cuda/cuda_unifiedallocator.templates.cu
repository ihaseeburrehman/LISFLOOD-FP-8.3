#include "cuda_unifiedallocator.cuh"
#include "cuda_util.cuh"

template<typename T>
T* lis::cuda::UnifiedAllocator<T>::allocate(std::size_t n)
{
	return static_cast<T*>(malloc_unified(sizeof(T) * n));
}

template<typename T>
void lis::cuda::UnifiedAllocator<T>::deallocate(T* ptr, std::size_t unused)
{
	free_unified(ptr);
}
