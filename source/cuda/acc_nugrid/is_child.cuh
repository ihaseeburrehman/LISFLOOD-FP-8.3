#pragma once

#include "cuda_runtime.h"
#include "MortonCode.h"
#include "index_1D.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__device__ __forceinline__ bool is_child
(
	MortonCode       fine_code,
	MortonCode       current_code,
	int              level,
	Parameters parameters
)

{
	return ( ( fine_code >> ( 2 * (parameters.L - level) ) ) == current_code);
}

}
}
}