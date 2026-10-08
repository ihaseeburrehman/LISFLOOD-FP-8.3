#pragma once
#include "../../lisflood.h"
#include "cuda_acc_flow.cuh"
#include "cuda_geometry.cuh"

namespace lis
{
namespace cuda
{
namespace acc
{

class MaxH
{
public:
	MaxH
	(
		Geometry& geometry
	);

	void update
	(
		const Flow& flow
	);

	~MaxH();

private:
	const int elements;
	void* temp;
	size_t bytes;
};

}
}
}
