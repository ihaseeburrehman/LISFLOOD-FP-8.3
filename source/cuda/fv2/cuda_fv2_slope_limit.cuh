#pragma once
#include "../../lisflood.h"
#include "cuda_fv2_flow.cuh"
#include "cuda_geometry.cuh"

namespace lis
{
namespace cuda
{
namespace fv2
{

struct Stencil
{
	NUMERIC_TYPE backward_const;
	NUMERIC_TYPE backward_slope;
	NUMERIC_TYPE local_const;
	NUMERIC_TYPE local_slope;
	NUMERIC_TYPE forward_const;
	NUMERIC_TYPE forward_slope;

	__device__
	NUMERIC_TYPE limit
	(
		NUMERIC_TYPE mesh_delta,
		NUMERIC_TYPE stencil_minH
	);

	__device__
	static NUMERIC_TYPE minH_x
	(
		Flow flow,
		int i,
		int j
	);

	__device__
	static NUMERIC_TYPE minH_y
	(
		Flow flow,
		int i,
		int j
	);

private:
	__device__ NUMERIC_TYPE backward_neg();
	__device__ NUMERIC_TYPE local_pos();
	__device__ NUMERIC_TYPE local_neg();
	__device__ NUMERIC_TYPE forward_pos();

};

struct Slopes
{
	NUMERIC_TYPE* __restrict__ ETA1x;
	NUMERIC_TYPE* __restrict__ ETA1y;
	NUMERIC_TYPE* __restrict__ HU1x;
	NUMERIC_TYPE* __restrict__ HU1y;
	NUMERIC_TYPE* __restrict__ HV1x;
	NUMERIC_TYPE* __restrict__ HV1y;

	__device__
	Stencil ETA_stencil_x
	(
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE* H,
		int i,
		int j
	);

	__device__
	Stencil HU_stencil_x
	(
		NUMERIC_TYPE* HU,
		int i,
		int j
	);

	__device__
	Stencil HV_stencil_x
	(
		NUMERIC_TYPE* HV,
		int i,
		int j
	);

	__device__
	Stencil ETA_stencil_y
	(
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE* H,
		int i,
		int j
	);

	__device__
	Stencil HU_stencil_y
	(
		NUMERIC_TYPE* HU,
		int i,
		int j
	);

	__device__
	Stencil HV_stencil_y
	(
		NUMERIC_TYPE* HV,
		int i,
		int j
	);

	static void allocate_device
	(
		Slopes& flow,
		Geometry& geometry
	);

	static void free_device
	(
		Slopes& slopes
	);
};

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
