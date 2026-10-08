#include "cuda_fv2_slope_limit.cuh"
#include "cuda_fv2_solver.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"
#include <helper_cuda.h>
#include <cub/cub.cuh>

__device__
NUMERIC_TYPE lis::cuda::fv2::Stencil::limit
(
	NUMERIC_TYPE mesh_delta,
	NUMERIC_TYPE stencil_minH
)
{
	if (stencil_minH <= C(0.09)) return C(0.0);

	NUMERIC_TYPE b_neg = backward_neg();
	NUMERIC_TYPE l_pos = local_pos();
	NUMERIC_TYPE l_neg = local_neg();
	NUMERIC_TYPE f_pos = forward_pos();

	NUMERIC_TYPE b_jump = l_pos - b_neg;
	NUMERIC_TYPE f_jump = l_neg - f_pos;

	NUMERIC_TYPE r;

	if (FABS(b_jump) <= C(0.0))
	{
		r = C(0.0);
	} 
	else 
	{
		r = f_jump / b_jump;
	}

	NUMERIC_TYPE beta = C(1.25);

	NUMERIC_TYPE phi = FMAX(C(0.0), FMAX(FMIN(beta * r, C(1.0)), FMIN(r, beta)));

	NUMERIC_TYPE limiter = phi * b_jump;

	return limiter;

}


__device__
NUMERIC_TYPE lis::cuda::fv2::Stencil::minH_x
(
	Flow flow,
	int i,
	int j
)
{
	NUMERIC_TYPE backward0 = flow.H[j*cuda::pitch + i-1];
	NUMERIC_TYPE local0 = flow.H[j*cuda::pitch + i];
	NUMERIC_TYPE forward0 = flow.H[j*cuda::pitch + i+1];

	return FMIN(backward0, FMIN(local0, forward0));
}

__device__
NUMERIC_TYPE lis::cuda::fv2::Stencil::minH_y
(
	Flow flow,
	int i,
	int j
)
{
	NUMERIC_TYPE backward0 = flow.H[(j+1)*cuda::pitch + i];
	NUMERIC_TYPE local0 = flow.H[j*cuda::pitch + i];
	NUMERIC_TYPE forward0 = flow.H[(j-1)*cuda::pitch + i];

	return FMIN(backward0, FMIN(local0, forward0));
}

__device__ NUMERIC_TYPE lis::cuda::fv2::Stencil::backward_neg()
{
	return backward_const;
}

__device__ NUMERIC_TYPE lis::cuda::fv2::Stencil::local_pos()
{
	return local_const;
}

__device__ NUMERIC_TYPE lis::cuda::fv2::Stencil::local_neg()
{
	return local_const;
}

__device__ NUMERIC_TYPE lis::cuda::fv2::Stencil::forward_pos()
{
	return forward_const;
}

__device__
lis::cuda::fv2::Stencil lis::cuda::fv2::Slopes::ETA_stencil_x
(
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* H,
	int i,
	int j
)
{
	return 
	{
		DEM[j * cuda::pitch + i - 1] + H[j * cuda::pitch + i - 1],
		0.0, //ETA1x[j * cuda::pitch + i - 1],
		DEM[j * cuda::pitch + i] + H[j * cuda::pitch + i],
		0.0, //ETA1x[j * cuda::pitch + i],
		DEM[j * cuda::pitch + i + 1] + H[j * cuda::pitch + i + 1],
		0.0 //ETA1x[j * cuda::pitch + i + 1]
	};
}

__device__
lis::cuda::fv2::Stencil lis::cuda::fv2::Slopes::HU_stencil_x
(
	NUMERIC_TYPE* HU,
	int i,
	int j
)
{
	return 
	{
		HU[j * cuda::pitch + i - 1],
		0.0, //HU1x[j * cuda::pitch + i - 1],
		HU[j * cuda::pitch + i],
		0.0, //HU1x[j * cuda::pitch + i],
		HU[j * cuda::pitch + i + 1],
		0.0, //HU1x[j * cuda::pitch + i + 1]
	};
}

__device__
lis::cuda::fv2::Stencil lis::cuda::fv2::Slopes::HV_stencil_x
(
	NUMERIC_TYPE* HV,
	int i,
	int j
)
{
	return 
	{
		HV[j * cuda::pitch + i - 1],
		0.0, //HV1x[j * cuda::pitch + i - 1],
		HV[j * cuda::pitch + i],
		0.0, //HV1x[j * cuda::pitch + i],
		HV[j * cuda::pitch + i + 1],
		0.0, //HV1x[j * cuda::pitch + i + 1]
	};
}

__device__
lis::cuda::fv2::Stencil lis::cuda::fv2::Slopes::ETA_stencil_y
(
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* H,
	int i,
	int j
)
{
	return 
	{
		DEM[(j + 1) * cuda::pitch + i] + H[(j + 1) * cuda::pitch + i],
		0.0, //ETA1y[(j + 1) * cuda::pitch + i],
		DEM[j * cuda::pitch + i] + H[j * cuda::pitch + i],
		0.0, //ETA1y[j * cuda::pitch + i],
		DEM[(j - 1) * cuda::pitch + i] + H[(j - 1) * cuda::pitch + i],
		0.0, //ETA1y[(j - 1) * cuda::pitch + i]
	};
}

__device__
lis::cuda::fv2::Stencil lis::cuda::fv2::Slopes::HU_stencil_y
(
	NUMERIC_TYPE* HU,
	int i,
	int j
)
{
	return 
	{
		HU[(j + 1) * cuda::pitch + i],
		0.0, //HU1y[(j + 1) * cuda::pitch + i],
		HU[j * cuda::pitch + i],
		0.0, //HU1y[j * cuda::pitch + i],
		HU[(j - 1) * cuda::pitch + i],
		0.0, //HU1y[(j - 1) * cuda::pitch + i]
	};
}

__device__
lis::cuda::fv2::Stencil lis::cuda::fv2::Slopes::HV_stencil_y
(
	NUMERIC_TYPE* HV,
	int i,
	int j
)
{
	return 
	{
		HV[(j + 1) * cuda::pitch + i],
		0.0, //HV1y[(j + 1) * cuda::pitch + i],
		HV[j * cuda::pitch + i],
		0.0, //HV1y[j * cuda::pitch + i],
		HV[(j - 1) * cuda::pitch + i],
		0.0, //HV1y[(j - 1) * cuda::pitch + i]
	};
}

void lis::cuda::fv2::Slopes::allocate_device
(
	Slopes& slopes,
	Geometry& geometry
)
{
	slopes.ETA1x = cuda::GhostRaster::allocate_device(geometry);
	slopes.ETA1y = cuda::GhostRaster::allocate_device(geometry);
	slopes.HU1x = cuda::GhostRaster::allocate_device(geometry);
	slopes.HU1y = cuda::GhostRaster::allocate_device(geometry);
	slopes.HV1x = cuda::GhostRaster::allocate_device(geometry);
	slopes.HV1y = cuda::GhostRaster::allocate_device(geometry);
}

void lis::cuda::fv2::Slopes::free_device
(
	Slopes& slopes
)
{
	cuda::free_device(slopes.ETA1x);
	cuda::free_device(slopes.ETA1y);
	cuda::free_device(slopes.HU1x);
	cuda::free_device(slopes.HU1y);
	cuda::free_device(slopes.HV1x);
	cuda::free_device(slopes.HV1y);
}

lis::cuda::fv2::MaxH::MaxH
(
	Geometry& geometry
)
:
elements(lis::GhostRaster::elements(geometry))
{
	temp = nullptr;
	NUMERIC_TYPE* dummy_in = nullptr;
	NUMERIC_TYPE* dummy_out = nullptr;

	checkCudaErrors(cub::DeviceReduce::Max(temp, bytes,
				dummy_in, dummy_out, elements));

	temp = malloc_device(bytes);
}

void lis::cuda::fv2::MaxH::update
(
	const Flow& flow
)
{
	checkCudaErrors(cub::DeviceReduce::Max(temp, bytes, flow.H,
				&(cuda::fv2::maxH), elements));
}

lis::cuda::fv2::MaxH::~MaxH()
{
	free_device(temp);
}
