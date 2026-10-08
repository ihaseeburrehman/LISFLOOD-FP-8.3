#include "../lisflood.h"
#include "cuda_unifiedallocator.cuh"
#include "cuda_unifiedallocator.templates.cu"

template class lis::cuda::UnifiedAllocator<NUMERIC_TYPE>;
