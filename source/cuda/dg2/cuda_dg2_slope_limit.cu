#include "cuda_dg2_slope_limit.cuh"
#include "cuda_dg2_solver.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"
#include <helper_cuda.h>
#include <cub/cub.cuh>

__device__
NUMERIC_TYPE lis::cuda::dg2::Stencil::limit
(
	NUMERIC_TYPE mesh_delta,
	NUMERIC_TYPE stencil_minH
)
{
	if (stencil_minH <= C(0.30)*cuda::dg2::maxH) return local_slope;

	NUMERIC_TYPE b_neg = backward_neg();
	NUMERIC_TYPE l_pos = local_pos();
	NUMERIC_TYPE l_neg = local_neg();
	NUMERIC_TYPE f_pos = forward_pos();

	NUMERIC_TYPE b_jump = FABS(l_pos - b_neg);
	NUMERIC_TYPE f_jump = FABS(l_neg - f_pos);

	NUMERIC_TYPE norm = FMAX(FABS(local_const-backward_const), FABS(forward_const-local_const));

	NUMERIC_TYPE DS_backward = C(0.0);
	NUMERIC_TYPE DS_forward = C(0.0);

	if (FABS(norm) > C(1e-6))
    {
        DS_backward = b_jump/(C(0.5)*mesh_delta*norm);
        DS_forward = f_jump/(C(0.5)*mesh_delta*norm);
    }

	if (FMAX(DS_backward, DS_forward) >= cuda::solver_params.KrivodonovaThresh)
	{
		return minmod
        (
            local_slope,
            (local_const-backward_const)/SQRT(C(3.0)),
            (forward_const-local_const)/SQRT(C(3.0))
        );
	}
	else
	{
		return local_slope;
	}
}

__device__
NUMERIC_TYPE lis::cuda::dg2::Stencil::minmod
(
	NUMERIC_TYPE a,
	NUMERIC_TYPE b,
	NUMERIC_TYPE c
)
{
    if (a*b > C(0.0) && a*c > C(0.0))
    {
        NUMERIC_TYPE sign_a = (a > C(0.0)) - (a < C(0.0));
        return sign_a * FMIN(FABS(a), FMIN(FABS(b), FABS(c)));
    }
    else
    {
        return C(0.0);
    }
}

__device__
NUMERIC_TYPE lis::cuda::dg2::Stencil::minH_x
(
	Flow flow,
	int i,
	int j
)
{
	NUMERIC_TYPE backward0 = flow.H[j*cuda::pitch + i-1];
	NUMERIC_TYPE backward1x = flow.H1x[j*cuda::pitch + i-1];
	NUMERIC_TYPE local0 = flow.H[j*cuda::pitch + i];
	NUMERIC_TYPE local1x = flow.H1x[j*cuda::pitch + i];
	NUMERIC_TYPE forward0 = flow.H[j*cuda::pitch + i+1];
	NUMERIC_TYPE forward1x = flow.H1x[j*cuda::pitch + i+1];

	NUMERIC_TYPE b_min = FMIN(backward0 - backward1x, backward0 + backward1x);
	NUMERIC_TYPE l_min = FMIN(local0 - local1x, local0 + local1x);
	NUMERIC_TYPE f_min = FMIN(forward0 - forward1x, forward0 + forward1x);

	return FMIN(b_min, FMIN(l_min, f_min));
}

__device__
NUMERIC_TYPE lis::cuda::dg2::Stencil::minH_y
(
	Flow flow,
	int i,
	int j
)
{
	NUMERIC_TYPE backward0 = flow.H[(j+1)*cuda::pitch + i];
	NUMERIC_TYPE backward1y = flow.H1y[(j+1)*cuda::pitch + i];
	NUMERIC_TYPE local0 = flow.H[j*cuda::pitch + i];
	NUMERIC_TYPE local1y = flow.H1y[j*cuda::pitch + i];
	NUMERIC_TYPE forward0 = flow.H[(j-1)*cuda::pitch + i];
	NUMERIC_TYPE forward1y = flow.H1y[(j-1)*cuda::pitch + i];

	NUMERIC_TYPE b_min = FMIN(backward0 - backward1y, backward0 + backward1y);
	NUMERIC_TYPE l_min = FMIN(local0 - local1y, local0 + local1y);
	NUMERIC_TYPE f_min = FMIN(forward0 - forward1y, forward0 + forward1y);

	return FMIN(b_min, FMIN(l_min, f_min));
}

__device__ NUMERIC_TYPE lis::cuda::dg2::Stencil::backward_neg()
{
	return backward_const + SQRT(C(3.0))*backward_slope;
}

__device__ NUMERIC_TYPE lis::cuda::dg2::Stencil::local_pos()
{
	return local_const - SQRT(C(3.0))*local_slope;
}

__device__ NUMERIC_TYPE lis::cuda::dg2::Stencil::local_neg()
{
	return local_const + SQRT(C(3.0))*local_slope;
}

__device__ NUMERIC_TYPE lis::cuda::dg2::Stencil::forward_pos()
{
	return forward_const - SQRT(C(3.0))*forward_slope;
}

__device__
lis::cuda::dg2::Stencil lis::cuda::dg2::Slopes::ETA_stencil_x
(
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* H,
	int i,
	int j
)
{
	return 
	{
		DEM[j*cuda::pitch + i-1] + H[j*cuda::pitch + i-1],
		ETA1x[j*cuda::pitch + i-1],
		DEM[j*cuda::pitch + i] + H[j*cuda::pitch + i],
		ETA1x[j*cuda::pitch + i],
		DEM[j*cuda::pitch + i+1] + H[j*cuda::pitch + i+1],
		ETA1x[j*cuda::pitch + i+1]
	};
}

__device__
lis::cuda::dg2::Stencil lis::cuda::dg2::Slopes::HU_stencil_x
(
	NUMERIC_TYPE* HU,
	int i,
	int j
)
{
	return 
	{
		HU[j*cuda::pitch + i-1],
		HU1x[j*cuda::pitch + i-1],
		HU[j*cuda::pitch + i],
		HU1x[j*cuda::pitch + i],
		HU[j*cuda::pitch + i+1],
		HU1x[j*cuda::pitch + i+1]
	};
}

__device__
lis::cuda::dg2::Stencil lis::cuda::dg2::Slopes::HV_stencil_x
(
	NUMERIC_TYPE* HV,
	int i,
	int j
)
{
	return 
	{
		HV[j*cuda::pitch + i-1],
		HV1x[j*cuda::pitch + i-1],
		HV[j*cuda::pitch + i],
		HV1x[j*cuda::pitch + i],
		HV[j*cuda::pitch + i+1],
		HV1x[j*cuda::pitch + i+1]
	};
}

__device__
lis::cuda::dg2::Stencil lis::cuda::dg2::Slopes::ETA_stencil_y
(
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* H,
	int i,
	int j
)
{
	return 
	{
		DEM[(j+1)*cuda::pitch + i] + H[(j+1)*cuda::pitch + i],
		ETA1y[(j+1)*cuda::pitch + i],
		DEM[j*cuda::pitch + i] + H[j*cuda::pitch + i],
		ETA1y[j*cuda::pitch + i],
		DEM[(j-1)*cuda::pitch + i] + H[(j-1)*cuda::pitch + i],
		ETA1y[(j-1)*cuda::pitch + i]
	};
}

__device__
lis::cuda::dg2::Stencil lis::cuda::dg2::Slopes::HU_stencil_y
(
	NUMERIC_TYPE* HU,
	int i,
	int j
)
{
	return 
	{
		HU[(j+1)*cuda::pitch + i],
		HU1y[(j+1)*cuda::pitch + i],
		HU[j*cuda::pitch + i],
		HU1y[j*cuda::pitch + i],
		HU[(j-1)*cuda::pitch + i],
		HU1y[(j-1)*cuda::pitch + i]
	};
}

__device__
lis::cuda::dg2::Stencil lis::cuda::dg2::Slopes::HV_stencil_y
(
	NUMERIC_TYPE* HV,
	int i,
	int j
)
{
	return 
	{
		HV[(j+1)*cuda::pitch + i],
		HV1y[(j+1)*cuda::pitch + i],
		HV[j*cuda::pitch + i],
		HV1y[j*cuda::pitch + i],
		HV[(j-1)*cuda::pitch + i],
		HV1y[(j-1)*cuda::pitch + i]
	};
}

void lis::cuda::dg2::Slopes::allocate_device
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

void lis::cuda::dg2::Slopes::free_device
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

lis::cuda::dg2::MaxH::MaxH
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

void lis::cuda::dg2::MaxH::update
(
	const Flow& flow
)
{
	checkCudaErrors(cub::DeviceReduce::Max(temp, bytes, flow.H,
				&(cuda::dg2::maxH), elements));
}

lis::cuda::dg2::MaxH::~MaxH()
{
	free_device(temp);
}
