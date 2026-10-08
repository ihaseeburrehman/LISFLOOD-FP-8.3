#pragma once
#include "../lisflood.h"
#include "cuda_flow.cuh"

namespace lis
{
namespace cuda
{

struct HLL
{
	__device__ static FlowVector x
	(
		const FlowVector& U_neg,
		const FlowVector& U_pos
	);

	__device__ static FlowVector y
	(
		const FlowVector& U_neg,
		const FlowVector& U_pos
	);
};

}
}
