#include "cuda_flow.cuh"
#include "cuda_solver.cuh"

__device__ lis::cuda::FlowVector lis::cuda::FlowVector::star
(
	NUMERIC_TYPE Z,
	NUMERIC_TYPE Zstar
) const
{
	NUMERIC_TYPE Hstar = calculate_Hstar(Z, Zstar);

	return { Hstar, Hstar*speed(HU), Hstar*speed(HV) };
};

__device__
lis::cuda::FlowVector lis::cuda::FlowVector::physical_flux_x() const
{
	if (H <= cuda::solver_params.DepthThresh)
	{
		return FlowVector();
	}
	else
	{
		return { HU, HU*HU/H + C(0.5)*cuda::physical_params.g*H*H, HU*HV/H };
	}
}

__device__
lis::cuda::FlowVector lis::cuda::FlowVector::physical_flux_y() const
{
	if (H <= cuda::solver_params.DepthThresh)
	{
		return FlowVector();
	}
	else
	{
		return { HV, HU*HV/H, HV*HV/H + C(0.5)*cuda::physical_params.g*H*H };
	}
}

__device__ NUMERIC_TYPE lis::cuda::FlowVector::calculate_Hstar
(
	NUMERIC_TYPE Z,
	NUMERIC_TYPE Zstar
) const
{
	NUMERIC_TYPE ETA = H + Z;
	return FMAX(C(0.0), ETA - Zstar);
}

__device__ NUMERIC_TYPE lis::cuda::FlowVector::speed
(
	NUMERIC_TYPE discharge
) const
{
	if (H >= cuda::solver_params.DepthThresh)
	{
		return discharge/H;
	}
	else
	{
		return C(0.0);
	}
}
