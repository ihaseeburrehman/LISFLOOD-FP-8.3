#include <helper_cuda.h>

template<typename T>
void lis::cuda::copy_to_symbol
(
	const T& symbol,
	const void* src,
	size_t count
)
{
	checkCudaErrors(cudaMemcpyToSymbol(symbol, src, count));
}

