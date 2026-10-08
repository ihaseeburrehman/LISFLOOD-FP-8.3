#include "cuda_hll.cuh"
#include "cuda_solver.cuh"

__device__ lis::cuda::FlowVector lis::cuda::HLL::y
(
	const FlowVector& u_neg,
	const FlowVector& u_pos
)
{
	FlowVector u_neg_rotated = { u_neg.H, u_neg.HV, -u_neg.HU };
	FlowVector u_pos_rotated = { u_pos.H, u_pos.HV, -u_pos.HU };

	FlowVector F_rotated = HLL::x(u_neg_rotated, u_pos_rotated);

	return { F_rotated.H, -F_rotated.HV, F_rotated.HU };
}

__device__ lis::cuda::FlowVector lis::cuda::HLL::x
(
	const FlowVector& u_neg,
	const FlowVector& u_pos
)
{
	if (u_neg.H <= cuda::solver_params.DepthThresh &&
			u_pos.H <= cuda::solver_params.DepthThresh)
	{
		return FlowVector();
	}

	NUMERIC_TYPE U_neg = C(0.0);
	NUMERIC_TYPE V_neg = C(0.0);
	if (u_neg.H > cuda::solver_params.DepthThresh)
	{
		U_neg = u_neg.HU / u_neg.H;
		V_neg = u_neg.HV / u_neg.H;
	}

	NUMERIC_TYPE U_pos = C(0.0); 
	NUMERIC_TYPE V_pos = C(0.0);
	if (u_pos.H > cuda::solver_params.DepthThresh)
	{
		U_pos = u_pos.HU / u_pos.H;
		V_pos = u_pos.HV / u_pos.H;
	}

	NUMERIC_TYPE A_neg = SQRT(cuda::physical_params.g * u_neg.H);
	NUMERIC_TYPE A_pos = SQRT(cuda::physical_params.g * u_pos.H);

	NUMERIC_TYPE H_star = POW(C(0.5)*(A_neg+A_pos) + C(0.25)*(U_neg - U_pos),
			C(2.0))	/ cuda::physical_params.g;
	NUMERIC_TYPE U_star = C(0.5)*(U_neg + U_pos) + A_neg - A_pos;
	NUMERIC_TYPE A_star = SQRT(cuda::physical_params.g*H_star);

	NUMERIC_TYPE S_neg;
	if (u_neg.H <= cuda::solver_params.DepthThresh)
	{
		S_neg = U_pos - C(2.0)*A_pos;
	}
	else
	{
		S_neg = FMIN(U_neg - A_neg, U_star - A_star);
	}

	NUMERIC_TYPE S_pos;
	if (u_pos.H <= cuda::solver_params.DepthThresh)
	{
		S_pos = U_neg + C(2.0)*A_neg;
	}
	else
	{
		S_pos = FMAX(U_pos + A_pos, U_star + A_star);
	}

	FlowVector F_neg = {
		u_neg.HU,
		U_neg * u_neg.HU + C(0.5) * cuda::physical_params.g * u_neg.H * u_neg.H,
		u_neg.HV*U_neg
	};

	FlowVector F_pos = {
		u_pos.HU,
		U_pos * u_pos.HU + C(0.5) * cuda::physical_params.g * u_pos.H * u_pos.H,
		u_pos.HV*U_pos
	};

	if (S_neg >= C(0.0))
	{
		return F_neg;
	}
	else if (S_neg < C(0.0) && S_pos >= C(0.0))
	{
		FlowVector F;

		F.H =
			(S_pos*F_neg.H - S_neg*F_pos.H + S_neg*S_pos*(u_pos.H-u_neg.H))
			/
			(S_pos - S_neg);
		F.HU =
			(S_pos*F_neg.HU - S_neg*F_pos.HU + S_neg*S_pos*(u_pos.HU-u_neg.HU))
			/
			(S_pos - S_neg);

		NUMERIC_TYPE S_mid =
			(S_neg*u_pos.H*(U_pos - S_pos) - S_pos*u_neg.H*(U_neg - S_neg))
			/
			(u_pos.H*(U_pos-S_pos) - u_neg.H*(U_neg - S_neg));

		if (S_neg < C(0.0) && S_mid >= C(0.0))
		{
			F.HV = F.H * V_neg;
		}
		else
		{
			F.HV = F.H * V_pos;
		}

		return F;
	}
	else
	{
		return F_pos;
	}
}
